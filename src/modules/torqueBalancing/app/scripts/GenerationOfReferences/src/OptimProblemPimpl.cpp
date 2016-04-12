#include "OptimProblemPimpl.h"

#include <yarp/sig/Vector.h>
#include <yarp/os/LogStream.h>
#include <IpTNLPAdapter.hpp>
#include <cassert>
#include <iDynTree/HighLevel/DynamicsComputations.h>
#include <iDynTree/Core/Transform.h>
#include <Eigen/Core>

using namespace Ipopt;

bool OptimProblemPimpl::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                                                   Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style)
{
    n = qDes.size();
    m = 3; // CoM position constraint
    if (feetInContact == BOTH_FEET_IN_CONTACT)
        m += 3 //relative transform position constraint
        + 4; //relative transform orientation (expressed in quaternion) constraint

    nnz_jac_g = m * n; //for now just treat everything as dense
    nnz_h_lag = n * n; //for now just treat everything as dense

    index_style = C_STYLE;
    return true;
}

bool OptimProblemPimpl::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                                                      Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u)
{
    for (Index i = 0; i < n; i++) {
        x_l[i] = -1e+19; //nlp_lower_bound_inf; //I cannot find the definition of this variable!
        x_u[i] =  1e+19;
    }
    //CoM constraints
    g_l[0] = g_u[0] = comDes[0];
    g_l[1] = g_u[1] = comDes[1];
    g_l[2] = g_u[2] = comDes[2];

    if (feetInContact == BOTH_FEET_IN_CONTACT) {
        // relative transform position constraint
        iDynTree::Position position = right_X_left.getPosition();
        iDynTree::Rotation rotation = right_X_left.getRotation();
        g_l[3] = g_u[3] = position(0);
        g_l[4] = g_u[4] = position(1);
        g_l[5] = g_u[5] = position(2);

        //TODO:        g_l[3] = g_u[3] = rotation.getQuaterion(0);
        g_l[6] = g_u[6] = rotation(0,0);
        g_l[7] = g_u[7] = rotation(0,0);
        g_l[8] = g_u[8] = rotation(0,0);
        g_l[9] = g_u[9] = rotation(0,0);
    }

    return true;
}

bool OptimProblemPimpl::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                                                         bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                                                         Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda)
{
    //Let's assert what we expect: only initial x_0
    assert(init_x);
    assert(!init_z);
    assert(!init_lambda);

    for (Index i = 0; i < n; ++i) {
        x[i] = qDes[i];
    }
    return true;
}

bool OptimProblemPimpl::eval_f(Ipopt::Index n, const Ipopt::Number* x,
                                             bool new_x, Ipopt::Number& obj_value)
{
    //What to do if new_x == false?
    //Computing objective: 1/2 * || q - q_des ||^2

    Eigen::Map<const Eigen::VectorXd> q(x, n);
    qError = q - qDes;
    obj_value = 0.5 * qError.transpose() * qError;
    return true;
}

bool OptimProblemPimpl::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                                                  Ipopt::Number* grad_f)
{
    //Gradient of the objective.
    //Simply q - qDes
    Eigen::Map<const Eigen::VectorXd> q(x, n);
    Eigen::Map<Eigen::VectorXd> gradient(grad_f, n);
    gradient = q - qDes;
    return true;
}

bool OptimProblemPimpl::eval_g(Ipopt::Index n, const Ipopt::Number* x,
                                             bool new_x, Ipopt::Index m, Ipopt::Number* g)
{
    for (unsigned index = 0; index < jointsMapping.size(); ++index) {
        allJoints(jointsMapping[index]) = qDes[index];
    }
    dynamics.setRobotState(allJoints, dofsSizeZero, dofsSizeZero, world_gravity);

    // CoM forward kinematic
    // must switch branch in iDynTree

    if (feetInContact == BOTH_FEET_IN_CONTACT) {

        iDynTree::Transform kinematic = dynamics.getRelativeTransform("l_sole", "r_sole");
        iDynTree::Position position = kinematic.getPosition();
        iDynTree::Rotation rotation = kinematic.getRotation();
        g[3] = position(0);
        g[4] = position(1);
        g[5] = position(2);

        //TODO:        g_l[3] = g_u[3] = rotation.getQuaterion(0);
        g[6] = rotation(0,0);
        g[7] = rotation(0,0);
        g[8] = rotation(0,0);
        g[9] = rotation(0,0);
    }
}

bool OptimProblemPimpl::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                                                 Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
                                                 Ipopt::Index *jCol, Ipopt::Number* values)
{
    if (!values) {
        //Sparsity structure of the Jacobian
    } else {
        //Actual Jacobian
        // CoM Jacobian
        // must switch branch in iDynTree

        if (feetInContact == BOTH_FEET_IN_CONTACT) {

        }
    }
}

bool OptimProblemPimpl::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                                             Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
                                             bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                                             Ipopt::Index* jCol, Ipopt::Number* values)
{
    //We can start with a quasi-Newton method to avoid computing the hessian of the constraints
    if (!values) {
        //Sparsity structure of the Hessian
    } else {
        //Actual Hessian

    }
}

void OptimProblemPimpl::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n,
                                                        const Ipopt::Number* x, const Ipopt::Number* z_L,
                                                        const Ipopt::Number* z_U, Ipopt::Index m, const Ipopt::Number* g,
                                                        const Ipopt::Number* lambda, Ipopt::Number obj_value,
                                                        const Ipopt::IpoptData* ip_data,
                                                        Ipopt::IpoptCalculatedQuantities* ip_cq)
{
}


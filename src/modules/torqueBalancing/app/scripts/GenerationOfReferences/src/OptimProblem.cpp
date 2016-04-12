#include "OptimProblem.h"
#include <iDynTree/HighLevel/DynamicsComputations.h>
#include <iDynTree/Core/Transform.h>
#include <yarp/sig/Vector.h>
#include <yarp/os/LogStream.h>
#include <IpTNLP.hpp>
#include <IpIpoptApplication.hpp>
#include <IpTNLPAdapter.hpp>
#include <cassert>
#include <Eigen/Core>

using namespace Ipopt;


#pragma mark - Private implementation

typedef enum _FeetInContact {
    LEFT_FOOT_IN_CONTACT,
    RIGHT_FOOT_IN_CONTACT,
    BOTH_FEET_IN_CONTACT,
} FeetInContact;

class OptimProblem::OptimProblemPimpl : public Ipopt::TNLP
{
public:
    iDynTree::HighLevel::DynamicsComputations dynamics;

    Eigen::VectorXd qDes;
    Eigen::VectorXd comDes;

    //handling feet in contact
    FeetInContact feetInContact;
    iDynTree::Transform leftToRightFootTransformation;

    //Buffers
    Eigen::VectorXd qError;
    Eigen::MatrixXd hessian;

#pragma mark - IpOpt methods

    virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                              Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style);

    virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                                 Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u);

    virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                                    bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                                    Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda);

    virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x,
                        bool new_x, Ipopt::Number& obj_value);

    virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                             Ipopt::Number* grad_f);

    virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x,
                        bool new_x, Ipopt::Index m, Ipopt::Number* g);

    virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                            Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
                            Ipopt::Index *jCol, Ipopt::Number* values);

    virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                        Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
                        bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                        Ipopt::Index* jCol, Ipopt::Number* values);

    virtual void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n,
                                   const Ipopt::Number* x, const Ipopt::Number* z_L,
                                   const Ipopt::Number* z_U, Ipopt::Index m, const Ipopt::Number* g,
                                   const Ipopt::Number* lambda, Ipopt::Number obj_value,
                                   const Ipopt::IpoptData* ip_data,
                                   Ipopt::IpoptCalculatedQuantities* ip_cq);


};

#pragma mark - OptimProblem implementation

OptimProblem::OptimProblem()
: pimpl(0)
{
    pimpl = new OptimProblemPimpl();
}

OptimProblem::~OptimProblem()
{
    if (pimpl) {
        delete pimpl;
        pimpl = 0;
    }
}

bool OptimProblem::initializeModel(std::string modelFile)
{
    return pimpl->dynamics.loadRobotModelFromFile(modelFile);
    //how to specify the joints we are interested in?
}


bool OptimProblem::solveOptimization(const yarp::sig::Vector& desiredCoM, const yarp::sig::Vector& desiredJoints, std::string feetInContact)
{
    Eigen::Map<const Eigen::VectorXd> inCom(desiredCoM.data(), desiredCoM.size());
    Eigen::Map<const Eigen::VectorXd> inQDes(desiredJoints.data(), desiredJoints.size());
    pimpl->qDes = inQDes;
    pimpl->comDes = inCom;
    if (feetInContact == "left") pimpl->feetInContact = LEFT_FOOT_IN_CONTACT;
    else if (feetInContact == "right") pimpl->feetInContact = RIGHT_FOOT_IN_CONTACT;
    else if (feetInContact == "both") {
        pimpl->feetInContact = BOTH_FEET_IN_CONTACT;
        //we must read the relative transform. I don't have the current state. I cannot compute it
        pimpl->leftToRightFootTransformation = pimpl->dynamics.getRelativeTransform("l_sole", "r_sole");
    } else {
        yError() << "Unsupported feet configuration";
        return false;
    }

    //resizing buffers
    pimpl->qError.resize(pimpl->qDes.size());
    pimpl->hessian.setIdentity(pimpl->qDes.size(), pimpl->qDes.size());

//    // Create a new instance of your nlp
//    //  (use a SmartPtr, not raw)
    SmartPtr<TNLP> mynlp = pimpl;
//
    // Create a new instance of IpoptApplication
    //  (use a SmartPtr, not raw)
    // We are using the factory, since this allows us to compile this
    // example with an Ipopt Windows DLL
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
//
//    // Change some options
//    // Note: The following choices are only examples, they might not be
//    //       suitable for your optimization problem.
//    app->Options()->SetNumericValue("tol", 1e-9);
//    app->Options()->SetStringValue("mu_strategy", "adaptive");
//    app->Options()->SetStringValue("output_file", "ipopt.out");
//
//    // Intialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Solve_Succeeded) {
        yError("*** Error during initialization of IpOpt!");
        return (int) status;
    }

    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(mynlp);

    if (status == Solve_Succeeded) {
        printf("\n\n*** The problem solved!\n");
    }
    else {
        printf("\n\n*** The problem FAILED!\n");
    }

    // As the SmartPtrs go out of scope, the reference count
    // will be decremented and the objects will automatically
    // be deleted.

    return (int) status;

}


#pragma mark - IpOpt methods

bool OptimProblem::OptimProblemPimpl::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
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

bool OptimProblem::OptimProblemPimpl::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
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
        iDynTree::Position position = leftToRightFootTransformation.getPosition();
        iDynTree::Rotation rotation = leftToRightFootTransformation.getRotation();
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

bool OptimProblem::OptimProblemPimpl::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
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

bool OptimProblem::OptimProblemPimpl::eval_f(Ipopt::Index n, const Ipopt::Number* x,
            bool new_x, Ipopt::Number& obj_value)
{
    //What to do if new_x == false?
    //Computing objective: 1/2 * || q - q_des ||^2

    Eigen::Map<const Eigen::VectorXd> q(x, n);
    qError = q - qDes;
    obj_value = 0.5 * qError.transpose() * qError;
    return true;
}

bool OptimProblem::OptimProblemPimpl::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                 Ipopt::Number* grad_f)
{
    //Gradient of the objective.
    //Simply q - qDes
    Eigen::Map<const Eigen::VectorXd> q(x, n);
    Eigen::Map<Eigen::VectorXd> gradient(grad_f, n);
    gradient = q - qDes;
    return true;
}

bool OptimProblem::OptimProblemPimpl::eval_g(Ipopt::Index n, const Ipopt::Number* x,
            bool new_x, Ipopt::Index m, Ipopt::Number* g)
{
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

bool OptimProblem::OptimProblemPimpl::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
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

bool OptimProblem::OptimProblemPimpl::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
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

void OptimProblem::OptimProblemPimpl::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n,
                               const Ipopt::Number* x, const Ipopt::Number* z_L,
                               const Ipopt::Number* z_U, Ipopt::Index m, const Ipopt::Number* g,
                               const Ipopt::Number* lambda, Ipopt::Number obj_value,
                               const Ipopt::IpoptData* ip_data,
                               Ipopt::IpoptCalculatedQuantities* ip_cq)
{
}


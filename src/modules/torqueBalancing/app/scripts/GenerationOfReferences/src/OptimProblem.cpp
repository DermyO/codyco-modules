#include "OptimProblem.h"
#include "OptimProblemPimpl.h"

#include <yarp/sig/Vector.h>
#include <yarp/os/LogStream.h>
#include <Eigen/Core>
#include <IpIpoptApplication.hpp>

using namespace Ipopt;

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

bool OptimProblem::initializeModel(const std::string modelFile, const std::vector<std::string> &jointsMapping)
{
    bool result = pimpl->dynamics.loadRobotModelFromFile(modelFile);
    if (!result) {
        yError() << "Error loading URDF model from " << modelFile;
        return false;
    }

    yInfo() << "Model loaded with " << pimpl->dynamics.getNrOfDegreesOfFreedom() << " DoFs";

    pimpl->allJoints.resize(pimpl->dynamics.getNrOfDegreesOfFreedom());
    pimpl->dofsSizeZero.resize(pimpl->dynamics.getNrOfDegreesOfFreedom());
    pimpl->allJoints.zero();
    pimpl->dofsSizeZero.zero();

    if (jointsMapping.empty())
        pimpl->jointsMapping.clear();
    else {
        pimpl->jointsMapping.reserve(jointsMapping.size());
        pimpl->qDes.setZero(jointsMapping.size());
        //mapping contains the list of joints (name) which are used in the optimization
        //Idea:
        // allJoints = initial (zero) values
        // before set state I use the content of jointsMapping: which is a map between the
        // optimization variable index (qDes) and the allJoints index.
        for (std::vector<std::string>::const_iterator it = jointsMapping.begin();
             it != jointsMapping.end(); ++it) {
            //TODO: getJointIndex missing @traversaro
            int index = pimpl->dynamics.getJointIndex(*it);
            pimpl->jointsMapping.push_back(index);
        }
        yInfo() << "Mapping translated to (# => model index): " << pimpl->jointsMapping;
    }

    //looking for feet frames
    pimpl->leftFootFrameID = pimpl->dynamics.getFrameIndex("l_sole");
    pimpl->rightFootFrameID = pimpl->dynamics.getFrameIndex("r_sole");
    return result;
}


bool OptimProblem::solveOptimization(const yarp::sig::Vector& desiredCoM, const yarp::sig::Vector& desiredJoints, std::string feetInContact)
{
    if (pimpl->jointsMapping.size() > 0)
        assert(desiredJoints.size() == pimpl->qDes.size());

    Eigen::Map<const Eigen::VectorXd> inCom(desiredCoM.data(), desiredCoM.size());
    Eigen::Map<const Eigen::VectorXd> inQDes(desiredJoints.data(), desiredJoints.size());
    pimpl->qDes = inQDes;
    pimpl->comDes = inCom;
    if (feetInContact == "left") pimpl->feetInContact = LEFT_FOOT_IN_CONTACT;
    else if (feetInContact == "right") pimpl->feetInContact = RIGHT_FOOT_IN_CONTACT;
    else if (feetInContact == "both") {
        pimpl->feetInContact = BOTH_FEET_IN_CONTACT;
        //we must read the relative transform. I don't have the current state. I cannot compute it
        pimpl->right_X_left = pimpl->dynamics.getRelativeTransform("l_sole", "r_sole");
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
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
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

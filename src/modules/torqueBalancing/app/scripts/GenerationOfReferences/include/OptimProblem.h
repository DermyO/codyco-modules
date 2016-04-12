#ifndef OPTIM_PROBLEM_H
#define OPTIM_PROBLEM_H

#include <IpTNLP.hpp>

namespace yarp {
    namespace sig {
        class Vector;
    }
}

class OptimProblem
{

    class OptimProblemPimpl;
    OptimProblemPimpl *pimpl;


public:

    OptimProblem();
    virtual ~OptimProblem();

    bool initializeModel(std::string modelFile);
    bool solveOptimization(const yarp::sig::Vector& desiredCoM, const yarp::sig::Vector& desiredJoints, std::string feetInContact);


};

#endif // OPTIM_PROBLEM_H

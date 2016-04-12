#ifndef OPTIM_PROBLEM_H
#define OPTIM_PROBLEM_H

#include <IpTNLP.hpp>
#include <vector>

namespace yarp {
    namespace sig {
        class Vector;
    }
}

class OptimProblemPimpl;

class OptimProblem
{
    OptimProblemPimpl *pimpl;


public:

    OptimProblem();
    virtual ~OptimProblem();

    bool initializeModel(const std::string modelFile, const std::vector<std::string>& jointsMapping = std::vector<std::string>());
    bool solveOptimization(const yarp::sig::Vector& desiredCoM, const yarp::sig::Vector& desiredJoints, std::string feetInContact);


};

#endif // OPTIM_PROBLEM_H

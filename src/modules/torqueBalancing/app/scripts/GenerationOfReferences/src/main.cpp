

#include <yarp/os/ResourceFinder.h>
#include <yarp/os/Value.h>
#include <yarp/os/LogStream.h>
#include <yarp/sig/Vector.h>
#include "OptimProblem.h"


/**
 *
 */
int main(int argc, char **argv) {
    using namespace yarp::os;

    ResourceFinder resourceFinder = ResourceFinder::getResourceFinderSingleton();
    resourceFinder.configure(argc, argv);

    if (resourceFinder.check("help")) {
        std::cout<< "Possible parameters" << std::endl << std::endl;
        std::cout<< "\t--context          :Where to find a user defined .ini file within $ICUB_ROOT/app e.g. /adaptiveControl/conf" << std::endl;
        std::cout<< "\t--robotURDFFile    :URDF file name" << std::endl;
        std::cout<< "\t--comDes           :Desired CoM of the robot." << std::endl;
        std::cout<< "\t--qDes             :Desired joint positions of the robot." << std::endl;
        std::cout<< "\t--feetInSupport    :left, right or both" << std::endl;
        return 0;
    }

    //read model file name
    std::string filename = resourceFinder.check("robotURDFFile", Value("model.urdf"), "Checking for model URDF file").asString();
    std::string filepath = resourceFinder.findFileByName(filename);

    yInfo() << "Robot model found in " << filepath;
    
    //read desired CoM
    if (!resourceFinder.check("comDes", "Checking desired CoM parameter")) {
        yError("Parameter comDes is required");
        return -1;
    }
    Value &comDes = resourceFinder.find("comDes");
    //Start checking: it should be a list of 3
    if (!comDes.isList() || comDes.asList()->size() != 3) {
        yError("Number of elements in comDes parameter is wrong. Expecting 3 values");
        return -2;
    }
    yarp::sig::Vector desiredCoM(3);
    Bottle* comList = comDes.asList();
    desiredCoM[0] = comList->get(0).asDouble();
    desiredCoM[1] = comList->get(1).asDouble();
    desiredCoM[2] = comList->get(2).asDouble();

    //read desired Joints configuration
    //(Which joints!? The model has more joints than we need. We have to map this stuff)
    if (!resourceFinder.check("qDes", "Checking desired joint configuration parameter")) {
        yError("Parameter qDes is required");
        return -1;
    }
    Value &qDes = resourceFinder.find("qDes");

    //read which foot/feet is in support
    if (!resourceFinder.check("feetInSupport", "Checking feet in support parameter")) {
        yError("Parameter feetInSupport is required");
        return -1;
    }
    Value &feet = resourceFinder.find("feetInSupport");


    OptimProblem problem;
    if (!problem.initializeModel(filepath)) {
        yError("Error initializing the robot model");
        return -2;
    }

    problem.solveOptimization(desiredCoM, desiredCoM, feet.asString());


    return 0;
}

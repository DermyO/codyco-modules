

#include <yarp/os/ResourceFinder.h>
#include <yarp/os/Value.h>
#include <yarp/os/LogStream.h>
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
    //read desired Joints configuration
    //(Which joints!? The model has more joints than we need. We have to map this stuff)
    //read which foot/feet is in support


    OptimProblem problem;
    problem.initializeModel(filepath);



    return 0;
}

/*
 * Copyright (C) 2014 CoDyCo
 * Author: Daniele Pucci, Francesco Romano
 * email:  daniele.pucci@iit.it, francesco.romano@iit.it
 * Permission is granted to copy, distribute, and/or modify this program
 * under the terms of the GNU General Public License, version 2 or any
 * later version published by the Free Software Foundation.
 *
 * A copy of the license can be found at
 * http://www.robotcub.org/icub/license/gpl.txt
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details
 */

/**
 *
 @ingroup codyco_module
 Copyright (C) 2014 CoDyCo Project
 
 CopyPolicy: Released under the terms of the GNU GPL v2.0.
 
 **/

/**
 Notes:
 There are two preprocessor macros: TORQUE_CONTROL  and GAZEBO_SIMULATOR.
 Torque control defines if I can use or not the ITorqueControl interface.
 If I have the torque control interface I can directly set the torque references to the robot (either the simulator or the real one)
 If I do not have the torque control interface (but this is a temporary situation... with t_available -> +inf) I implement the low level (raw) torque control, by setting directly the voltages to the motors.
 
 Currently (but it is late) I don't know if the GAZEBO_SIMULATOR is needed anymore....
 
 */

#ifndef ADAPTIVECONTROLTHREAD_H
#define ADAPTIVECONTROLTHREAD_H

#include "AdaptiveControlConstants.h"
#include <yarp/os/RateThread.h>
#include <yarp/sig/Vector.h>
#include <string>

#include <Eigen/Core>

namespace paramHelp {
    class ParamHelperServer;
    class ParamHelperClient;
}

namespace yarp {
    namespace dev {
        class PolyDriver;
        class IEncodersTimed;
        class IControlMode;
#ifdef ADAPTIVECONTROL_TORQUECONTROL
        class ITorqueControl;
#else
//        class ITorqueControlRaw;
#endif
    }
    namespace os {
        class Bottle;
        template <class T> class BufferedPort;
    }
}

namespace iCub {
    namespace ctrl {
        class AWLinEstimator;
    }
}

namespace adaptiveControl {
    
	const int passiveJointIndex = 0;
	const int activeJointIndex = 3;
	
#ifndef ADAPTIVECONTROL_TORQUECONTROL
	//class MotorParameters;
#endif
	
    class AdaptiveControlThread: public yarp::os::RateThread, public paramHelp::CommandObserver {
        
    private:
        //internal state variables
        bool _controlEnabled;
        unsigned short _maxReadFailed;
        unsigned short _failedReads;
        bool _firstRunLoop;
        double _initialTime;
        double _previousTime;
        bool _piHatModificationOn;
        
        
        //configuration parameters
        const std::string &_threadName;
        const std::string &_robotName;
        const std::string &_robotPart;
        paramHelp::ParamHelperServer &_paramServer;
        
        //in-out varables
        yarp::dev::PolyDriver* _driver;
        yarp::dev::IEncodersTimed* _encoders;
        yarp::dev::IControlMode* _controlMode;
#ifdef ADAPTIVECONTROL_TORQUECONTROL
        yarp::dev::ITorqueControl* _torqueControl;
#else
//        yarp::dev::ITorqueControlRaw* _rawTorqueControl;
#endif
        yarp::os::BufferedPort<yarp::os::Bottle>* _torqueOutput;
		yarp::os::BufferedPort<yarp::sig::Vector>* _debugPort;
        
        //Temp: to read directly from gazebo
        yarp::os::BufferedPort<yarp::os::Bottle>* _speedInput;
        
#ifndef ADAPTIVECONTROL_TORQUECONTROL
        paramHelp::ParamHelperClient& _paramClient; //used to send torques commands to torque control
#endif
		
        iCub::ctrl::AWLinEstimator* _velocityEstimator;
        int _outputEnabled;
        
        
        double _minDeterminantValue;
        
        //reference trajectory: for now i compute it internally. In the future we can read from a port
        double _refBaseline;
        double _refAngularVelocity; //in rad/s
        double _refAmplitude;
        double _refPhase; //in rad
        
        double _currentRef;
        
        //geometric parameters
        double _link1Length;
        double _link2Length;
        
        //gains
        double _lambda;
        Eigen::Vector2d _kappa;
        Eigen::Matrix8d _Gamma;
        
        //Variables for computing control
        Eigen::Vector2d _q;
        Eigen::Vector2d _dq;
        Eigen::Vector2d _xi;
        Eigen::Vector8d _piHat;
        
        //variables update rules
        Eigen::Vector2d _dxi;
        Eigen::Vector8d _dpiHat;
        
        //Streaming output parameters
        yarp::sig::Vector _outputTau;
        
        //debug variables
        double _massMatrixDeterminant;
   
        yarp::dev::PolyDriver* openDriver(std::string localName, std::string robotName, std::string bodyPartName);
        void computeRegressor(const Eigen::Vector2d& q, /* Joint positions*/
                                          const Eigen::Vector2d& dq, /* Joint velocities*/
                                          const Eigen::Vector2d& dq_lin, /* Joint velocities. This is the term which multiplies linearly the C(.) term */
                                          const Eigen::Vector2d& ddq, /* Joint accelerations*/
                              Eigen::Matrix28d& regressor); /* output variable */
        bool readSensors(Eigen::Vector2d& positions, Eigen::Vector2d& velocities);
        void computeControl();
        void writeOutputs();
        void startControl();
        void stopControl();
		
		void writeDebug();
        
        void dampedPseudoInverse28(const Eigen::Matrix28d &A,
                                 double dampingFactor,
                                 Eigen::Matrix<double, 8, 2> &Apinv);
         void dampedPseudoInverse88(const Eigen::Matrix<double, 8, 8> &A,
                                 double dampingFactor,
                                 Eigen::Matrix<double, 8, 8> &Apinv);
		
    public:
        AdaptiveControlThread(const std::string& threadName,
                              const std::string& robotName,
                              const std::string& robotPart,
                              int periodMilliseconds,
                              paramHelp::ParamHelperServer&paramHelperServer,
#ifndef ADAPTIVECONTROL_TORQUECONTROL
                              paramHelp::ParamHelperClient& paramHelperClient,
#endif
                              const Eigen::Vector2d &linklengths);
        ~AdaptiveControlThread();
        
        bool setInitialConditions(const Eigen::Vector8d& initialPiHat, const double& initialXi1);
        bool controlEnabled();
        
        /* Overrided functions (Rate Thread)*/
        bool threadInit();
        void run();
        void threadRelease();
        
        void commandReceived(const paramHelp::CommandDescription &cd, const yarp::os::Bottle &params, yarp::os::Bottle &reply);
        
    };
}

#endif

/*
 * Copyright (C) 2016 Fondazione Istituto Italiano di Tecnologia - Italian Institute of Technology
 * Author: Jorhabib Eljaik
 * email:  jorhabib.eljaik@iit.it
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

#ifndef ICUB_WALKING_IK_THREAD_H
#define ICUB_WALKING_IK_THREAD_H

#include <yarp/os/RateThread.h>
#include <yarp/os/Property.h>
#include <yarp/os/BufferedPort.h>
#include <yarp/sig/Vector.h>
#include <yarp/os/ResourceFinder.h>
#include <yarp/os/Time.h>

#include <yarpWholeBodyInterface/yarpWholeBodyModel.h>
#include <yarpWholeBodyInterface/yarpWholeBodyStates.h>

#include "SplineInterpolator.h"
#include "UtilityFunctions.h"
#include "IK.h"

struct walkingParams{
    double z_c;
    int    n_strides;
    int    T_stride;
    int    T_switch;
    double step_width;
    double step_length;
    double step_height;
    int    n_samples;
    double g;
};

class iCubWalkingIKThread: public yarp::os::RateThread {
private:
    int m_period;
    std::string m_walkingPatternFile;
    walkingParams m_walkingParams;
    yarp::os::ResourceFinder m_rf;
    wbi::iWholeBodyModel* m_wbm;
    wbi::iWholeBodyStates* m_wbs;
    std::string m_outputDir;

public:
    iCubWalkingIKThread (int period,
                         wbi::iWholeBodyModel* wbm,
                         wbi::iWholeBodyStates* wbs,
                         walkingParams params,
                         yarp::os::ResourceFinder& rf,
                         std::string walkingPatternFile,
                         std::string outputDir);
    bool threadInit();
    void run();
    void threadRelease();
    void generateFeetTrajectories(std::string walkingPatternFile,
                                  walkingParams params);
    void inverseKinematics(walkingParams params);
};

#endif
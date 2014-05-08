/**
 * Copyright (C) 2014 CoDyCo
 * @author: Francesco Romano
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

#include "TorqueBalancingController.h"
#include "Reference.h"
#include <wbi/wholeBodyInterface.h>
#include <wbi/wbiUtil.h>
#include <codyco/MathUtils.h>
#include <codyco/LockGuard.h>


#include <Eigen/LU>
#include <Eigen/SVD>

namespace codyco {
    namespace torquebalancing {
        
        TorqueBalancingController::TorqueBalancingController(int period, ControllerReferences& references, wbi::wholeBodyInterface& robot)
        : RateThread(period)
        , m_robot(robot)
        , m_active(false)
        , m_centerOfMassLinkID(wbi::wholeBodyInterface::COM_LINK_ID)
        , m_references(references)
        , m_desiredJointsConfiguration(actuatedDOFs)
        , m_internal_desiredJointsConfiguration(actuatedDOFs)
        , m_centroidalMomentumGain(0)
        , m_internal_centroidalMomentumGain(0)
        , m_impedanceGains(actuatedDOFs)
        , m_internal_impedanceGains(actuatedDOFs)
        , m_desiredCOMAcceleration(3)
        , m_desiredFeetForces(12)
        , m_desiredCentroidalMomentum(6)
        , m_jointPositions(totalDOFs)
        , m_jointVelocities(totalDOFs)
        , m_torques(actuatedDOFs)
        , m_baseVelocity(6)
        , m_centerOfMassPosition(3)
        , m_rightFootPosition(7)
        , m_leftFootPosition(7)
        , m_feetJacobian(6 + 6, totalDOFs)
        , m_feetDJacobianDq(12)
        , m_massMatrix(totalDOFs, totalDOFs)
        , m_generalizedBiasForces(totalDOFs)
        , m_gravityBiasTorques(totalDOFs)
        , m_centroidalMomentum(6)
        , m_pseudoInverseOfJcMInvSt(actuatedDOFs, 12)
        , m_centroidalForceMatrix(6, 12)
        , m_gravityForce(6)
        , m_torquesSelector(totalDOFs, actuatedDOFs)
        , m_rotoTranslationVector(7) {}
        
        TorqueBalancingController::~TorqueBalancingController() {}

#pragma mark - RateThread methods
        bool TorqueBalancingController::threadInit()
        {
            using namespace Eigen;
            //Initialize constant variables
            m_robot.getLinkId("l_sole", m_leftFootLinkID);
            m_robot.getLinkId("r_sole", m_rightFootLinkID);
            
            m_leftFootToBaseRotationFrame.R = wbi::Rotation(0, 0, 1,
                                                            0, -1, 0,
                                                            1, 0, 0);
            
            //centroidal force matrix
            m_centroidalForceMatrix.setZero();
            m_centroidalForceMatrix.block<3, 3>(0, 0) = Matrix3d::Identity();
            m_centroidalForceMatrix.block<3, 3>(0, 6) = Matrix3d::Identity();
            m_centroidalForceMatrix.block<3, 3>(3, 3) = Matrix3d::Identity();
            m_centroidalForceMatrix.block<3, 3>(3, 9) = Matrix3d::Identity();
            //gravity
            m_gravityForce.setZero();
            m_gravityUnitVector[0] = m_gravityUnitVector[1] = 0;
            //TODO: check sign of gravity
            m_gravityUnitVector[2] = -9.81;
            return true;
        }
        
        void TorqueBalancingController::threadRelease()
        {
            
        }
        
        void TorqueBalancingController::run()
        {
            if (!this->isActiveState()) return;
            //read references
            readReferences();
            
            //read / update state
            updateRobotState();
            
            //compute desired feet forces
            computeFeetForces(m_desiredCOMAcceleration, m_desiredFeetForces);
            
            //compute torques
            computeTorques(m_desiredFeetForces, m_torques);
            
            //write torques
            writeTorques();
            
        }
        
#pragma mark - Getter and setter
        
        double TorqueBalancingController::centroidalMomentumGain()
        {
            codyco::LockGuard guard(m_mutex);
            return m_centroidalMomentumGain;
        }
        
        void TorqueBalancingController::setCentroidalMomentumGain(double centroidalMomentumGain)
        {
            codyco::LockGuard guard(m_mutex);
            m_centroidalMomentumGain = centroidalMomentumGain;
        }
        
        const Eigen::VectorXd& TorqueBalancingController::impedanceGains()
        {
            codyco::LockGuard guard(m_mutex);
            return m_impedanceGains;
        }
        
        void TorqueBalancingController::setImpedanceGains(Eigen::VectorXd& impedanceGains)
        {
            codyco::LockGuard guard(m_mutex);
            m_impedanceGains = impedanceGains;
        }
        
        const Eigen::VectorXd& TorqueBalancingController::desiredJointsConfiguration()
        {
            codyco::LockGuard guard(m_mutex);
            return m_desiredJointsConfiguration;
        }
        
        void TorqueBalancingController::setDesiredJointsConfiguration(Eigen::VectorXd& desiredJointsConfiguration)
        {
            codyco::LockGuard guard(m_mutex);
            m_desiredJointsConfiguration = desiredJointsConfiguration;
        }
        
        void TorqueBalancingController::setActiveState(bool isActive)
        {
            codyco::LockGuard guard(m_mutex);
            if (m_active == isActive) return;
            if (isActive) {
                m_desiredCOMAcceleration.setZero(); //reset reference
                m_robot.setControlMode(wbi::CTRL_MODE_TORQUE);
            } else {
                m_robot.setControlMode(wbi::CTRL_MODE_POS);
            }
            m_active = isActive;
        }
        
        bool TorqueBalancingController::isActiveState()
        {
            codyco::LockGuard guard(m_mutex);
            return m_active;
        }
        
#pragma mark - Controller methods
        
        void TorqueBalancingController::readReferences()
        {
            if (m_references.desiredCOMAcceleration().isValid())
                m_desiredCOMAcceleration = m_references.desiredCOMAcceleration().value();
            codyco::LockGuard guard(m_mutex);
            m_internal_centroidalMomentumGain = m_centroidalMomentumGain;
            m_internal_desiredJointsConfiguration = m_internal_desiredJointsConfiguration;
            m_internal_impedanceGains = m_impedanceGains;
        }
        
        bool TorqueBalancingController::updateRobotState()
        {
            //update world to base frame
            m_robot.computeH(m_jointPositions.data(), wbi::Frame(), m_leftFootLinkID, m_world2BaseFrame);
            m_world2BaseFrame = m_world2BaseFrame * m_leftFootToBaseRotationFrame;
            m_world2BaseFrame.setToInverse();
            
            //read positions and velocities
            m_robot.getEstimates(wbi::ESTIMATE_JOINT_POS, m_jointPositions.data());
            m_robot.getEstimates(wbi::ESTIMATE_JOINT_VEL, m_jointVelocities.data());
            
            //update mass matrix
            m_robot.computeMassMatrix(m_jointPositions.data(), m_world2BaseFrame, m_massMatrix.data());
            m_robot.computeCentroidalMomentum(m_jointPositions.data(), m_world2BaseFrame, m_jointVelocities.data(), m_baseVelocity.data(), m_centroidalMomentum.data());

            m_robot.forwardKinematics(m_jointPositions.data(), m_world2BaseFrame, m_centerOfMassLinkID, m_rotoTranslationVector.data());
            m_centerOfMassPosition = m_rotoTranslationVector.head<3>();
            m_robot.forwardKinematics(m_jointPositions.data(), m_world2BaseFrame, m_leftFootLinkID, m_leftFootPosition.data());
            m_robot.forwardKinematics(m_jointPositions.data(), m_world2BaseFrame, m_rightFootLinkID, m_rightFootPosition.data());
            
            //update jacobians (both feet in one variable)
            m_robot.computeJacobian(m_jointPositions.data(), m_world2BaseFrame, m_leftFootLinkID, m_feetJacobian.topRows(6).data());
            m_robot.computeJacobian(m_jointPositions.data(), m_world2BaseFrame, m_rightFootLinkID, m_feetJacobian.bottomRows(6).data());
            
            m_robot.computeDJdq(m_jointPositions.data(), m_world2BaseFrame, m_jointVelocities.data(), m_baseVelocity.data(), m_leftFootLinkID, m_feetDJacobianDq.head(6).data());
            m_robot.computeDJdq(m_jointPositions.data(), m_world2BaseFrame, m_jointVelocities.data(), m_baseVelocity.data(), m_rightFootLinkID, m_feetDJacobianDq.tail(6).data());

            
            //Compute bias forces
            m_robot.computeGeneralizedBiasForces(m_jointPositions.data(), m_world2BaseFrame, m_jointVelocities.data(), m_baseVelocity.data(), m_gravityUnitVector, m_generalizedBiasForces.data());
            m_robot.computeGeneralizedBiasForces(m_jointPositions.data(), m_world2BaseFrame, 0, 0, m_gravityUnitVector, m_gravityBiasTorques.data());

            
            return true;
        }
        
        void TorqueBalancingController::computeFeetForces(const Eigen::Ref<Eigen::MatrixXd>& desiredCOMAcceleration, Eigen::Ref<Eigen::MatrixXd> desiredFeetForces)
        {
            using namespace Eigen;
            double mass = m_massMatrix(0, 0);
            m_gravityForce(2) = -mass * 9.81;
            
            //building centroidalForceMatrix
            skewSymmentricMatrix(m_leftFootPosition.head<3>() - m_centerOfMassPosition, m_centroidalForceMatrix.block<3, 3>(3, 0));
            skewSymmentricMatrix(m_rightFootPosition.head<3>() - m_centerOfMassPosition, m_centroidalForceMatrix.block<3, 3>(3, 6));
            
            m_desiredCentroidalMomentum.head<3>() = mass * desiredCOMAcceleration;
            m_desiredCentroidalMomentum.tail<3>() = -m_internal_centroidalMomentumGain * m_centroidalMomentum.tail<3>();

            desiredFeetForces = m_centroidalForceMatrix.jacobiSvd(ComputeThinU | ComputeThinV).solve(m_desiredCentroidalMomentum - m_gravityForce);
            //TODO: we can also test LDLT decomposition since it requires positive semidefinite matrix
//            desiredFeetForces = m_centroidalForceMatrix.ldlt().solve(m_desiredCentroidalMomentum - m_gravityForce);
            
        }
        
        void TorqueBalancingController::computeTorques(const Eigen::Ref<Eigen::MatrixXd>& desiredFeetForces, Eigen::Ref<Eigen::MatrixXd> torques)
        {
            using namespace Eigen;
            
            //TODO: decide later if there is a performance benefit in moving the declaration of the variables in the class
            //Names are taken from "math" from brevity
            MatrixXd JcMInv = m_feetJacobian * m_massMatrix.inverse();
            MatrixXd JcMInvTorqueSelector = JcMInv * m_torquesSelector;
            
            math::pseudoInverse(JcMInvTorqueSelector, PseudoInverseTolerance, m_pseudoInverseOfJcMInvSt);
            MatrixXd nullSpaceProjector = MatrixXd::Identity(actuatedDOFs, actuatedDOFs) - m_pseudoInverseOfJcMInvSt * JcMInvTorqueSelector;
                        
            m_torques = m_pseudoInverseOfJcMInvSt * (JcMInv * m_generalizedBiasForces - m_feetDJacobianDq - JcMInv * m_feetJacobian.transpose() * desiredFeetForces);
            
            VectorXd torques0 = m_gravityBiasTorques.tail(actuatedDOFs) - m_feetJacobian.block(0, 6, totalDOFs, actuatedDOFs).transpose() * desiredFeetForces - m_internal_impedanceGains.asDiagonal() * (m_jointPositions.tail(actuatedDOFs) - m_internal_desiredJointsConfiguration);
            
            m_torques += nullSpaceProjector * torques0;
            
        }
        
        void TorqueBalancingController::writeTorques()
        {
            m_robot.setControlReference(m_torques.data());
        }
        
#pragma mark - Auxiliary functions
        
        void TorqueBalancingController::skewSymmentricMatrix(const Eigen::Ref<const Eigen::Vector3d>& vector, Eigen::Ref<Eigen::Matrix3d> skewSymmetricMatrix)
        {
            skewSymmetricMatrix.setZero();
//            S = [   0,   -w(3),    w(2);
//                 w(3),   0,     -w(1);
//                 -w(2),  w(1),     0   ];
            skewSymmetricMatrix(0, 1) = -vector(2);
            skewSymmetricMatrix(0, 2) = vector(1);
            skewSymmetricMatrix(1, 2) = -vector(0);
            skewSymmetricMatrix.bottomLeftCorner<2, 2>() = -skewSymmetricMatrix.topRightCorner<2, 2>().transpose();
        }
        
    }
}
/**********************************************************************************************
 * SixDOF - Version 0.0.1
 * by Henrik Söderlund <henrik.a.soderlund@gmail.com>
 *
 * Copyright (C) 2018 Henrik Söderlund
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 **********************************************************************************************/


#ifndef SixDOF_h
#define SixDOF_h

#include "MatrixMath.h"

#define INVERSE_JACOBIAN_METHOD
//#define PSEUDO_INVERSE_JACOBIAN_METHOD
//#define TRANSPOSE_JACOBIAN_METHOD

#define SINGULARITY_THRESHOLD 0.00001
#define REACHED_DESTINATION_THRESHOLD 0.005
#define INVKIN_TIMEOUT 1000

class SixDOF
{
  public:

    enum IKState {
        FAILED = -1,
        SUCCESS = 0,
        DEST_REACHED = 2
    };
    typedef enum IKState IKState_t;

    SixDOF(double* thetas, double* ds, double* alphas, double* as, int len);
    
    void forwardKinematics(double* angles);
    IKState_t inverseKinematics(double x, double y, double z, double phi, double theta, double psi, double grippingOffset, double alpha);
    IKState_t inverseKinematics2(double x, double y, double z, double phi, double theta, double psi, double grippingOffset, double alpha);

    void getPose(double* returnPose);
    void getJointAngles(double* returnAngles);

    IKState_t getIKStatus();

  private:
  
    int numJoints = 6; // Defaults to 6
    int numLinks = 7; // Defaults to 7

    //Manipulator pose
    mtx_type pose[6];

    // Unit vectors
    mtx_type v_i[3];
    mtx_type v_j[3];
    mtx_type v_k[3];

    // Identity matrix
    mtx_type eye[4][4];
    mtx_type eye3[3][3];

    IKState_t ikStatus;
 
    // Joint struct (Revolute)
    typedef struct 
    {
        // Position (global) and angle of joint
        double x;
        double y;
        double z;
        double a;
        
        // DH Parameters of joint
        double DH_Params[4];
        
        // Joint Jacobian
        mtx_type J[6]; // The size does not correspond to numJoints!
        // Homogenous transform matrix from fromJoint to toJoint
        mtx_type H[4][4]; 
    } Joint;
    
    typedef struct 
    {
        double length;
        
        Joint* fromJoint;
        Joint* toJoint;
        
        bool isEndEffector = false;
        
        // Homogenous transform matrix from origin to fromJoint
        mtx_type H0[4][4]; 
    } Link;
    
    // List of links that makes up the manipulator
    Link* links;

    IKState_t inverseWristKinematics(double x, double y, double z, int wristCenter, double alpha);
    
    void computeDH(double theta, double d, double alpha, double a, mtx_type* H);
    void computeJointJacobian(mtx_type* endH, mtx_type* prevH, mtx_type* J);
    void screwZ(double theta, double d, mtx_type* H);
    void screwX(double alpha, double a, mtx_type* H);
    void crossProduct(mtx_type* u, mtx_type* v, mtx_type* w);
    double angleDiff(mtx_type* v, mtx_type* u);
    double vecDistance(mtx_type* from, mtx_type* to);
    double norm(mtx_type* v);
    //double norm6(mtx_type* v);
    double dot(mtx_type* u, mtx_type* v);
    double determinant(mtx_type* M, int len);

    void rot2axis(mtx_type* R, mtx_type* ax);
    //void axis2rot(mtx_type* ax, mtx_type* R);
    //void eigDecomp(mtx_type* A, int m, int n, int iterations, mtx_type* out);
    
    void rotZ(double theta, mtx_type* R);
    void rotY(double phi, mtx_type* R);
    void rotX(double psi, mtx_type* R);
   
    void rot2euler(mtx_type* R, mtx_type* euler);
    void euler2rot(mtx_type* euler, mtx_type* R);
    void t2r(mtx_type* H, mtx_type* R);
    double traceR(mtx_type* R);
    double maxDiagR(mtx_type* R);
    int maxDiagRCol(mtx_type* R);
    //void euler2angvel(mtx_type* eulerFrom, mtx_type* eulerTo, mtx_type* angVel);
    
};

#endif

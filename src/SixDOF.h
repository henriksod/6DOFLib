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
#define REACHED_DESTINATION_THRESHOLD 0.01

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
    //IKState_t inverseKinematics(double x, double y, double z, double ax, double ay, double az, double theta, double alpha);

    IKState_t inverseKinematics(double x, double y, double z, double phi, double theta, double psi, double alpha);
    //IKState_t inverseKinematics(double x, double y, double z, double phi, double theta, double psi);
    
    void getPose(double* returnPose);
    void getJointAngles(double* returnAngles);

    IKState_t getIKStatus();

  private:
  
    int numJoints = 6; // Defaults to 6
    int numLinks = 7; // Defaults to 7
    
    mtx_type bufferV[3];
    mtx_type bufferU[3];
    mtx_type bufferW[3];
    mtx_type bufferX[3];
    mtx_type bufferY[3];
    mtx_type bufferZ[3];
    mtx_type bufferH[4][4];
    mtx_type bufferH2[4][4];
    mtx_type bufferR[3][3];
    mtx_type bufferR2[3][3];
    mtx_type bufferR3[3][3];
    mtx_type bufferR4[3][3];
    mtx_type bufferR5[3][3];
    
    mtx_type bufferAX[4];
    mtx_type bufferAY[4];

    mtx_type oc[3];

    mtx_type eye[4][4];
    
    // Manipulator Jacobian
    mtx_type* J;
    mtx_type* invJ;
    mtx_type* bufferinvJ;
    mtx_type* transposeJ;
    mtx_type* bufferJ;

    mtx_type bufferPose[6];
    
    //Manipulator pose
    mtx_type pose[6];
    mtx_type newPose[6];
    
    //Joint angle differences
    mtx_type* dAngles;

    // Unit vectors
    mtx_type v_i[3];
    mtx_type v_j[3];
    mtx_type v_k[3];

    mtx_type zPrev[3];
    mtx_type oPrev[3];
    mtx_type oEnd[3];

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
    
    void computeDH(double theta, double d, double alpha, double a, mtx_type* H);
    void computeJointJacobian(mtx_type* endH, mtx_type* prevH, mtx_type* J);
    void screwZ(double theta, double d, mtx_type* H);
    void screwX(double alpha, double a, mtx_type* H);
    void crossProduct(mtx_type* u, mtx_type* v, mtx_type* w);
    double angleDiff(mtx_type* v, mtx_type* u);
    double vecDistance(mtx_type* from, mtx_type* to);
    double determinant(mtx_type* M, int len);

    /*
    void rot2axis(mtx_type* R, mtx_type* ax);
    void axis2rot(mtx_type* ax, mtx_type* R);
    void eigDecomp(mtx_type* A, int m, int n, int iterations, mtx_type* out);
    */

    
    void rotZ(double theta, mtx_type* R);
    void rotY(double phi, mtx_type* R);
    void rotX(double psi, mtx_type* R);
   
    void rot2euler(mtx_type* R, mtx_type* euler);
    void euler2rot(mtx_type* euler, mtx_type* R);
    void euler2angvel(mtx_type* eulerFrom, mtx_type* eulerTo, mtx_type* angVel);
    
};

#endif

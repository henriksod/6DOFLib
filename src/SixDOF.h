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

#define JACOBIAN_SINGULARITY_THRESHOLD 0.00001
#define ALPHA 0.1

class SixDOF
{
  public:

    SixDOF(float* thetas, float* ds, float* alphas, float* as, int len);
    
    void forwardKinematics(float* angles);
    int inverseKinematics(float x, float y, float z, float theta, float phi, float psi, float dt);
    
    void getPose(float* returnPose);
    void getJointAngles(float* returnAngles);

  private:
  
    int numJoints = 6; // Defaults to 6
    int numLinks = 7; // Defaults to 7
    
    mtx_type bufferV[3];
    mtx_type bufferU[3];
    mtx_type bufferH[4][4];
    mtx_type bufferH2[4][4];

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
 
    // Joint struct (Revolute)
    typedef struct 
    {
        // Position (global) and angle of joint
        float x;
        float y;
        float z;
        float a;
        
        // DH Parameters of joint
        float DH_Params[4];
        
        // Joint Jacobian
        mtx_type J[6]; // The size does not correspond to numJoints!
        // Homogenous transform matrix from fromJoint to toJoint
        mtx_type H[4][4]; 
    } Joint;
    
    typedef struct 
    {
        float length;
        
        Joint* fromJoint;
        Joint* toJoint;
        
        bool isEndEffector = false;
        
        // Homogenous transform matrix from origin to fromJoint
        mtx_type H0[4][4]; 
    } Link;
    
    // List of links that makes up the manipulator
    Link* links;
    
    void computeDH(float theta, float d, float alpha, float a, mtx_type* H);
    void computeJointJacobian(mtx_type* endH, mtx_type* prevH, mtx_type* J);
    void screwZ(float theta, float d, mtx_type* H);
    void screwX(float alpha, float a, mtx_type* H);
    void crossProduct(mtx_type* u, mtx_type* v, mtx_type* w);
    float angleDiff(mtx_type* v, mtx_type* u);
    float determinant(mtx_type* M, int len);
    
};

#endif

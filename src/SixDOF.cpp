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


#include <stdlib.h>
#include "Arduino.h"
#include "SixDOF.h"


SixDOF::SixDOF(double* thetas, double* ds, double* alphas, double* as, int len)
{
    this->numJoints = len;
    this->numLinks = this->numJoints+1;
    
    this->links = (Link*)malloc(sizeof(Link)*(this->numLinks));

    for (int i = 0; i < this->numLinks-1; i++)
    {
        this->links[i].fromJoint = (Joint*)malloc(sizeof(Joint));
        
        this->links[i].fromJoint->x = 0;
        this->links[i].fromJoint->y = 0;
        this->links[i].fromJoint->z = 0;
        this->links[i].fromJoint->a = 0;
        
        this->links[i].fromJoint->DH_Params[0] = thetas[i];
        this->links[i].fromJoint->DH_Params[1] = ds[i];
        this->links[i].fromJoint->DH_Params[2] = alphas[i];
        this->links[i].fromJoint->DH_Params[3] = as[i];
        
        this->links[i].length = ds[i] == 0 ? as[i] : ds[i];
        this->links[i].isEndEffector = false;

        if (i > 0) {
            this->links[i-1].toJoint = (Joint*)malloc(sizeof(Joint));
            this->links[i-1].toJoint = links[i].fromJoint;
        }
    }

    this->links[this->numLinks-2].fromJoint = links[this->numLinks-3].toJoint;
    this->links[this->numLinks-2].isEndEffector = true;
    
    this->eye[0][0] = 1; this->eye[0][1] = 0; this->eye[0][2] = 0; this->eye[0][3] = 0;
    this->eye[1][0] = 0; this->eye[1][1] = 1; this->eye[1][2] = 0; this->eye[1][3] = 0;
    this->eye[2][0] = 0; this->eye[2][1] = 0; this->eye[2][2] = 1; this->eye[2][3] = 0;
    this->eye[3][0] = 0; this->eye[3][1] = 0; this->eye[3][2] = 0; this->eye[3][3] = 1;

    this->J          = (mtx_type*)malloc(6 * this->numJoints * sizeof(mtx_type));
    this->invJ       = (mtx_type*)malloc(this->numJoints * 6 * sizeof(mtx_type));
    this->bufferJ    = (mtx_type*)malloc(6 * this->numJoints * sizeof(mtx_type));
    this->bufferinvJ = (mtx_type*)malloc(this->numJoints * this->numJoints * sizeof(mtx_type));
    this->transposeJ = (mtx_type*)malloc(this->numJoints * 6 * sizeof(mtx_type));
    
    this->dAngles = (mtx_type*)malloc(sizeof(mtx_type)*this->numJoints);
    
    this->pose[0] = 0;
    this->pose[1] = 0;
    this->pose[2] = 0;
    this->pose[3] = 0;
    this->pose[4] = 0;
    this->pose[5] = 0;

    this->v_i[0] = 1;
    this->v_i[1] = 0;
    this->v_i[2] = 0;
    this->v_j[0] = 0;
    this->v_j[1] = 1;
    this->v_j[2] = 0;
    this->v_k[0] = 0;
    this->v_k[1] = 0;
    this->v_k[2] = 1;

    this->ikStatus = SUCCESS;
}

void SixDOF::getPose(double* returnPose)
{
    Link endEffector = this->links[this->numJoints-1];
    
    this->pose[0] = endEffector.H0[0][3];
    this->pose[1] = endEffector.H0[1][3];
    this->pose[2] = endEffector.H0[2][3];

    
    returnPose[0] = (double)this->pose[0];
    returnPose[1] = (double)this->pose[1];
    returnPose[2] = (double)this->pose[2];
    returnPose[3] = (double)this->pose[3];
    returnPose[4] = (double)this->pose[4];
    returnPose[5] = (double)this->pose[5];
}

void SixDOF::getJointAngles(double* returnAngles)
{
    for (int i = 1; i < this->numLinks; i++)
    {
        returnAngles[i-1] = this->links[i-1].fromJoint->a;
    }
}

void SixDOF::forwardKinematics(double* angles)
{
    Matrix.Copy((mtx_type*)this->eye, 4, 4, (mtx_type*)this->links[0].H0);
    Matrix.Copy((mtx_type*)this->eye, 4, 4, (mtx_type*)this->links[0].fromJoint->H);
    double theta = this->links[0].fromJoint->DH_Params[0];
    double d     = this->links[0].fromJoint->DH_Params[1];
    double alpha = this->links[0].fromJoint->DH_Params[2];
    double a     = this->links[0].fromJoint->DH_Params[3];
    computeDH(angles[0]+theta, d, alpha, a, (mtx_type*)this->links[0].H0);
    computeDH(angles[0]+theta, d, alpha, a, (mtx_type*)this->links[0].fromJoint->H);
    //Matrix.Print((mtx_type*)this->links[0].H0, 4, 4, "H0");

    this->links[0].fromJoint->x = this->links[0].H0[0][3];
    this->links[0].fromJoint->y = this->links[0].H0[1][3];
    this->links[0].fromJoint->z = this->links[0].H0[2][3];
    
    this->links[0].fromJoint->a = angles[0];
    
    for (int i = 1; i < this->numJoints; i++)
    {
        Matrix.Copy((mtx_type*)this->links[i-1].H0, 4, 4, (mtx_type*)this->links[i].H0);
        Matrix.Copy((mtx_type*)this->eye, 4, 4, (mtx_type*)this->links[i].fromJoint->H);
        
        theta = this->links[i].fromJoint->DH_Params[0];
        d     = this->links[i].fromJoint->DH_Params[1];
        alpha = this->links[i].fromJoint->DH_Params[2];
        a     = this->links[i].fromJoint->DH_Params[3];
        computeDH(angles[i]+theta, d, alpha, a, (mtx_type*)this->links[i].H0);
        computeDH(angles[i]+theta, d, alpha, a, (mtx_type*)this->links[i].fromJoint->H);
        
        this->links[i].fromJoint->x = this->links[i].H0[0][3];
        this->links[i].fromJoint->y = this->links[i].H0[1][3];
        this->links[i].fromJoint->z = this->links[i].H0[2][3];
        
        this->links[i].fromJoint->a = angles[i];
    }

    Link endEffector = this->links[this->numLinks-2];
    
    this->pose[0] = endEffector.H0[0][3];
    this->pose[1] = endEffector.H0[1][3];
    this->pose[2] = endEffector.H0[2][3];

    this->bufferR[0][0] = endEffector.H0[0][0]; this->bufferR[0][1] = endEffector.H0[0][1]; this->bufferR[0][2] = endEffector.H0[0][2];
    this->bufferR[1][0] = endEffector.H0[1][0]; this->bufferR[1][1] = endEffector.H0[1][1]; this->bufferR[1][2] = endEffector.H0[1][2];
    this->bufferR[2][0] = endEffector.H0[2][0]; this->bufferR[2][1] = endEffector.H0[2][1]; this->bufferR[2][2] = endEffector.H0[2][2];
    
    rot2euler((mtx_type*)bufferR, (mtx_type*)bufferV);
    this->pose[3] = bufferV[0];
    this->pose[4] = bufferV[1];
    this->pose[5] = bufferV[2];
}

////////////////////////////////////////// KINEMATIC DECOUPLING, IS IT POSSIBLE FOR N DOF? //////////////////////////////////////////
/*
SixDOF::IKState_t SixDOF::inverseKinematics(double x, double y, double z, double phi, double theta, double psi, double alpha)
{
    int wristCenter = this->numJoints-2;
    for (int i = this->numJoints-2; i > 0; i--)
    {
        if (this->links[i].length != 0) {
          wristCenter = i; 
          break;
        }
    }

    this->bufferV[0] = x; 
    this->bufferV[1] = y; 
    this->bufferV[2] = z; 

    this->bufferU[0] = 0; 
    this->bufferU[1] = 0; 
    this->bufferU[2] = 1; 

    rotX(phi, (mtx_type*)this->bufferR3);
    rotY(theta, (mtx_type*)this->bufferR3);
    rotZ(psi, (mtx_type*)this->bufferR3);

    Matrix.Multiply((mtx_type*)this->bufferR3, (mtx_type*)this->bufferU, 3, 3, 1, (mtx_type*)this->bufferW);
    Matrix.Scale((mtx_type*)this->bufferW, 3, 1, this->links[this->numJoints-1].length);
    Matrix.Subtract((mtx_type*)this->bufferV, (mtx_type*)this->bufferW, 3, 1, (mtx_type*)this->oc);

    this->bufferV[0] = this->links[wristCenter].fromJoint->x;
    this->bufferV[1] = this->links[wristCenter].fromJoint->y;
    this->bufferV[2] = this->links[wristCenter].fromJoint->z;

    double baseRotation = atan2(y,x);

    int numSolve = wristCenter-1;
    mtx_type* Jv = (mtx_type*)malloc(3 * numSolve * sizeof(mtx_type));
    mtx_type* JvT = (mtx_type*)malloc(3 * numSolve * sizeof(mtx_type));
    mtx_type* invJv = (mtx_type*)malloc(numSolve * numSolve * sizeof(mtx_type));
    mtx_type* pinvJv = (mtx_type*)malloc(3 * numSolve * sizeof(mtx_type));

    Link wristCenterLink = this->links[numSolve];

    //// WHILE...
    for (int i = 1; i < wristCenter; i++)
    {
        numSolve++;
        computeJointJacobian((mtx_type*)wristCenterLink.H0, (mtx_type*)this->links[i-1].H0, (mtx_type*)this->links[i].fromJoint->J);
        Jv[numSolve * 0 + i] = this->links[i].fromJoint->J[0];
        Jv[numSolve * 1 + i] = this->links[i].fromJoint->J[1];
        Jv[numSolve * 2 + i] = this->links[i].fromJoint->J[2];
    }

    Matrix.Transpose((mtx_type*)Jv, 3, numSolve, (mtx_type*)JvT);
    Matrix.Multiply((mtx_type*)JvT, (mtx_type*)Jv, numSolve, 3, numSolve, (mtx_type*)invJv);
    
    int success = Matrix.Invert((mtx_type*)invJv, numSolve);
    if (success == 0) {
      free(Jv);
      free(JvT);
      free(invJv);
      free(pinvJv);
      this->ikStatus = FAILED;
      return FAILED;
    }

    Matrix.Multiply((mtx_type*)invJv, (mtx_type*)JvT, numSolve, numSolve, 3, (mtx_type*)pinvJv);
    
    double SEDistance = vecDistance((mtx_type*)this->bufferV, (mtx_type*)this->oc);
    Matrix.Subtract((mtx_type*)this->oc, (mtx_type*)this->bufferV, 3, 1, (mtx_type*)this->bufferU);

    mtx_type* diffAngles = (mtx_type*)malloc(numSolve * sizeof(mtx_type));
    Matrix.Multiply((mtx_type*)pinvJv, (mtx_type*)this->bufferU, numSolve, 3, 1, (mtx_type*)diffAngles);

    mtx_type* newAngles = (mtx_type*)malloc((numSolve+1) * sizeof(mtx_type));
    newAngles[0] = baseRotation;
    for (int i = 1; i < numSolve; i++)
    {
        newAngles[i] = this->links[i+1].fromJoint->a + alpha*diffAngles[i];
    }
    
    forwardKinematics((double*)newAngles, numSolve+1);

    wristCenterLink = this->links[numSolve];
    this->bufferR[0][0] = wristCenterLink.H0[0][0]; this->bufferR[0][1] = wristCenterLink.H0[0][1]; this->bufferR[0][2] = wristCenterLink.H0[0][2];
    this->bufferR[1][0] = wristCenterLink.H0[1][0]; this->bufferR[1][1] = wristCenterLink.H0[1][1]; this->bufferR[1][2] = wristCenterLink.H0[1][2];
    this->bufferR[2][0] = wristCenterLink.H0[2][0]; this->bufferR[2][1] = wristCenterLink.H0[2][1]; this->bufferR[2][2] = wristCenterLink.H0[2][2];

    if (SEDistance <= REACHED_DESTINATION_THRESHOLD) {
      this->ikStatus = DEST_REACHED;
      return DEST_REACHED; 
    }

    free(diffAngles);
    free(newAngles);
    free(Jv);
    free(JvT);
    free(invJv);
    free(pinvJv);
    this->ikStatus = SUCCESS;
    return SUCCESS;
}
*/

//////////////////////////////////////////


SixDOF::IKState_t SixDOF::inverseKinematics(double x, double y, double z, double phi, double theta, double psi, double alpha)
{
    Link endEffector = this->links[this->numLinks-2];
    
    this->ikStatus = SUCCESS;
    
    this->newPose[0] = x;
    this->newPose[1] = y;
    this->newPose[2] = z;
    this->newPose[3] = phi;
    this->newPose[4] = theta;
    this->newPose[5] = psi;

    for(int i = 0; i < 6; i++)
          for(int j = 0; j < this->numJoints; j++)
                this->J[this->numJoints*i + j] = 0;

    
    computeJointJacobian((mtx_type*)endEffector.H0, (mtx_type*)this->eye, (mtx_type*)this->links[0].fromJoint->J);
    ////Matrix.Print((mtx_type*)this->eye, 4, 4, "eye");
    //Matrix.Print((mtx_type*)endEffector.H0, 4, 4, "H0");
    //Matrix.Print((mtx_type*)this->links[0].fromJoint->J, 6, 1, "J0");
    this->J[this->numJoints * 0 + 0] = this->links[0].fromJoint->J[0];
    this->J[this->numJoints * 1 + 0] = this->links[0].fromJoint->J[1];
    this->J[this->numJoints * 2 + 0] = this->links[0].fromJoint->J[2];
    this->J[this->numJoints * 3 + 0] = this->links[0].fromJoint->J[3];
    this->J[this->numJoints * 4 + 0] = this->links[0].fromJoint->J[4];
    this->J[this->numJoints * 5 + 0] = this->links[0].fromJoint->J[5];
    
    for (int i = 1; i < this->numLinks-1; i++)
    {
        computeJointJacobian((mtx_type*)endEffector.H0, (mtx_type*)this->links[i-1].H0, (mtx_type*)this->links[i].fromJoint->J);
        //Matrix.Print((mtx_type*)this->links[i-1].H0, 4, 4, "H");
        //Matrix.Print((mtx_type*)this->links[i].fromJoint->J, 6, 1, "J"); 
        
        this->J[this->numJoints * 0 + i] = this->links[i].fromJoint->J[0];
        this->J[this->numJoints * 1 + i] = this->links[i].fromJoint->J[1];
        this->J[this->numJoints * 2 + i] = this->links[i].fromJoint->J[2];
        this->J[this->numJoints * 3 + i] = this->links[i].fromJoint->J[3];
        this->J[this->numJoints * 4 + i] = this->links[i].fromJoint->J[4];
        this->J[this->numJoints * 5 + i] = this->links[i].fromJoint->J[5];
    }
    
    // Compute determinant of J > SINGULARITY_THRESHOLD
    double detJ = determinant((mtx_type*)this->J, this->numJoints);
    if (abs(detJ) > SINGULARITY_THRESHOLD)
    {
        int success = 0;
#ifdef INVERSE_JACOBIAN_METHOD
        if (this->numJoints == 6) {
          Matrix.Copy((mtx_type*)this->J, this->numJoints, this->numJoints, (mtx_type*)this->invJ);
          success = Matrix.Invert((mtx_type*)this->invJ, this->numJoints);
          if (success == 0) {
            this->ikStatus = FAILED;
            return FAILED;
          }
        } else {
          Matrix.Transpose((mtx_type*)this->J, 6, this->numJoints, (mtx_type*)this->transposeJ);
          Matrix.Multiply((mtx_type*)this->transposeJ, (mtx_type*)this->J, this->numJoints, 6, this->numJoints, (mtx_type*)this->bufferinvJ);
          
          success = Matrix.Invert((mtx_type*)this->bufferinvJ, this->numJoints);
          if (success == 0) {
            this->ikStatus = FAILED;
            return FAILED;
          }

          Matrix.Multiply((mtx_type*)this->bufferinvJ, (mtx_type*)this->transposeJ, this->numJoints, this->numJoints, 6, (mtx_type*)this->invJ);
        }
#endif

#ifdef PSEUDO_INVERSE_JACOBIAN_METHOD
        Matrix.Transpose((mtx_type*)this->J, 6, this->numJoints, (mtx_type*)this->transposeJ);
        Matrix.Multiply((mtx_type*)this->transposeJ, (mtx_type*)this->J, this->numJoints, 6, this->numJoints, (mtx_type*)this->bufferinvJ);
        
        success = Matrix.Invert((mtx_type*)this->bufferinvJ, this->numJoints);
        if (success == 0) {
          this->ikStatus = FAILED;
          return FAILED;
        }

        Matrix.Multiply((mtx_type*)this->bufferinvJ, (mtx_type*)this->transposeJ, this->numJoints, this->numJoints, 6, (mtx_type*)this->invJ);
#endif

#ifdef TRANSPOSE_JACOBIAN_METHOD
        Matrix.Transpose((mtx_type*)this->J, 6, this->numJoints, (mtx_type*)this->invJ);
#endif

        // SE
        this->bufferV[0] = this->pose[0];
        this->bufferV[1] = this->pose[1];
        this->bufferV[2] = this->pose[2];
        this->bufferU[0] = this->newPose[0];
        this->bufferU[1] = this->newPose[1];
        this->bufferU[2] = this->newPose[2];
        double SEDistance = vecDistance((mtx_type*)bufferV, (mtx_type*)bufferU);
        
        Matrix.Subtract((mtx_type*)this->bufferU, (mtx_type*)this->bufferV, 3, 1, (mtx_type*)this->bufferW);
        this->bufferPose[0] = this->bufferW[0];
        this->bufferPose[1] = this->bufferW[1];
        this->bufferPose[2] = this->bufferW[2];

        // SO
        this->bufferX[0] = this->pose[3];
        this->bufferX[1] = this->pose[4];
        this->bufferX[2] = this->pose[5];
        this->bufferY[0] = this->newPose[3];
        this->bufferY[1] = this->newPose[4];
        this->bufferY[2] = this->newPose[5];

        euler2angvel((mtx_type*)this->bufferX, (mtx_type*)this->bufferY, (mtx_type*)this->bufferZ);

        this->bufferPose[3] = this->bufferZ[0];
        this->bufferPose[4] = this->bufferZ[1];
        this->bufferPose[5] = this->bufferZ[2];
        
        Matrix.Multiply((mtx_type*)this->invJ, (mtx_type*)this->bufferPose, this->numJoints, 6, 1, (mtx_type*)this->dAngles);

        for (int i = 0; i < this->numJoints; i++)
        {
            this->dAngles[i] = this->links[i].fromJoint->a + alpha*this->dAngles[i];
        }
        
        forwardKinematics((double*)this->dAngles);

        if (SEDistance <= REACHED_DESTINATION_THRESHOLD) {
          this->ikStatus = DEST_REACHED;
          return DEST_REACHED; 
        }
    } else { 
      Serial.println("Jacobian singular - Determinant below defined JACOBIAN_SINGULARITY_THRESHOLD");
      this->ikStatus = FAILED;
      return FAILED; 
    }

    this->ikStatus = SUCCESS;
    return SUCCESS;
}

SixDOF::IKState_t SixDOF::getIKStatus()
{
  return this->ikStatus;
}

void SixDOF::computeDH(double theta, double d, double alpha, double a, mtx_type* H) 
{
    screwZ(theta, d, H);
    screwX(alpha, a, H);
}

void SixDOF::computeJointJacobian(mtx_type* endH, mtx_type* prevH, mtx_type* J) 
{
    // m rows, n cols
    
    this->zPrev[0] = prevH[4 * 0 + 2];
    this->zPrev[1] = prevH[4 * 1 + 2];
    this->zPrev[2] = prevH[4 * 2 + 2];
    
    this->oPrev[0] = prevH[4 * 0 + 3];
    this->oPrev[1] = prevH[4 * 1 + 3];
    this->oPrev[2] = prevH[4 * 2 + 3];
    
    this->oEnd[0] = endH[4 * 0 + 3];
    this->oEnd[1] = endH[4 * 1 + 3];
    this->oEnd[2] = endH[4 * 2 + 3];

    Matrix.Subtract((mtx_type*)this->oEnd, (mtx_type*)this->oPrev, 3, 1, (mtx_type*)this->bufferU);
    crossProduct(this->zPrev, this->bufferU, this->bufferV);
    J[0] = this->bufferV[0];
    J[1] = this->bufferV[1];
    J[2] = this->bufferV[2];
    J[3] = this->zPrev[0];
    J[4] = this->zPrev[1];
    J[5] = this->zPrev[2];
}

void SixDOF::screwZ(double theta, double d, mtx_type* H)
{
    double c=cos(theta);
    double s=sin(theta);

    this->bufferH[0][0] = c;  this->bufferH[0][1] = -s; this->bufferH[0][2] = 0; this->bufferH[0][3] = 0;
    this->bufferH[1][0] = s; this->bufferH[1][1] = c; this->bufferH[1][2] = 0; this->bufferH[1][3] = 0;
    this->bufferH[2][0] = 0;  this->bufferH[2][1] = 0; this->bufferH[2][2] = 1; this->bufferH[2][3] = d;
    this->bufferH[3][0] = 0;  this->bufferH[3][1] = 0; this->bufferH[3][2] = 0; this->bufferH[3][3] = 1;
    
    Matrix.Multiply((mtx_type*)H, (mtx_type*)this->bufferH, 4, 4, 4, (mtx_type*)this->bufferH2);
    Matrix.Copy((mtx_type*)this->bufferH2, 4, 4, (mtx_type*)H);
}

void SixDOF::screwX(double alpha, double a, mtx_type* H)
{
    double c=cos(alpha);
    double s=sin(alpha);

    this->bufferH[0][0] = 1; this->bufferH[0][1] = 0; this->bufferH[0][2] = 0;  this->bufferH[0][3] = a;
    this->bufferH[1][0] = 0; this->bufferH[1][1] = c; this->bufferH[1][2] = -s; this->bufferH[1][3] = 0;
    this->bufferH[2][0] = 0; this->bufferH[2][1] = s; this->bufferH[2][2] = c;  this->bufferH[2][3] = 0;
    this->bufferH[3][0] = 0; this->bufferH[3][1] = 0; this->bufferH[3][2] = 0;  this->bufferH[3][3] = 1;
    
    Matrix.Multiply((mtx_type*)H, (mtx_type*)this->bufferH, 4, 4, 4, (mtx_type*)this->bufferH2);
    Matrix.Copy((mtx_type*)this->bufferH2, 4, 4, (mtx_type*)H);
}

void SixDOF::crossProduct(mtx_type* v, mtx_type* u, mtx_type* w) 
{
    double jk = v[1]*u[2]-u[1]*v[2];
    double ik = v[0]*u[2]-u[0]*v[2];
    double ij = v[0]*u[1]-u[0]*v[1];
    
    Matrix.Scale((mtx_type*)this->v_i, 3, 1, jk);
    Matrix.Scale((mtx_type*)this->v_j, 3, 1, ik);
    Matrix.Scale((mtx_type*)this->v_k, 3, 1, ij);
    
    Matrix.Subtract((mtx_type*)this->v_i, (mtx_type*)this->v_j, 3, 1, (mtx_type*)this->bufferV);
    Matrix.Add((mtx_type*)this->bufferV, (mtx_type*)this->v_k, 3, 1, (mtx_type*)w);

    this->v_i[0] = 1;
    this->v_i[1] = 0;
    this->v_i[2] = 0;
    this->v_j[0] = 0;
    this->v_j[1] = 1;
    this->v_j[2] = 0;
    this->v_k[0] = 0;
    this->v_k[1] = 0;
    this->v_k[2] = 1;
}


double SixDOF::angleDiff(mtx_type* v, mtx_type* u)
{
    double dotProd = v[0]*u[0] + v[1]*u[1] + v[2]*u[2];
    double normV = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    double normU = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    
    double c_ang = dotProd/(normV*normU);
    double s_ang = sqrt(1-c_ang*c_ang);
    
    return atan2(s_ang, c_ang);
} 


double SixDOF::vecDistance(mtx_type* from, mtx_type* to)
{
  double diffX = to[0]-from[0];
  double diffY = to[1]-from[1];
  double diffZ = to[2]-from[2];
  return sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ);
}

double SixDOF::determinant(mtx_type* M, int len) 
{
  double det = 0;
  if (len == 2) 
  {
      det = (double)(M[len * 0 + 0]*M[len * 1 + 1]-M[len * 0 + 1]*M[len * 1 + 0]);
  }
  else if (len > 2) 
  {
      int splitSize = (int)sqrt(len*len-2*len+1);
      mtx_type* temp_M = (mtx_type*)malloc(splitSize * splitSize * sizeof(mtx_type)); 
      
      for (int i = 0; i < len; i++)
      {
          int a = 0;
          for (int j = 1; j < len; j++) {
              int b = 0;
              for (int k = 0; k < len; k++) {
                  if (k != i) {
                      temp_M[splitSize * a + b] = M[len * j + k];
                      b++;
                  }
              }
              a++;
          }
          
          int o = -1;
          for (int j = 0; j <= i; j++) o = o*o;
          
          if (i == 0) {
              det = (double)(M[len * 0 + i]*o*determinant((mtx_type*)temp_M, splitSize));
          } else {
              det = det + (double)(M[len * 0 + i]*o*determinant((mtx_type*)temp_M, splitSize));
          }
      }
      free((mtx_type*)temp_M);
  }
  return det;
}

void SixDOF::rotZ(double theta, mtx_type* R)
{
    double c=cos(theta);
    double s=sin(theta);

    this->bufferR[0][0] = c;  this->bufferR[0][1] = -s; this->bufferR[0][2] = 0;
    this->bufferR[1][0] = s; this->bufferR[1][1] = c; this->bufferR[1][2] = 0;
    this->bufferR[2][0] = 0;  this->bufferR[2][1] = 0; this->bufferR[2][2] = 1;
    
    Matrix.Multiply((mtx_type*)R, (mtx_type*)this->bufferR, 3, 3, 3, (mtx_type*)this->bufferR2);
    Matrix.Copy((mtx_type*)this->bufferR2, 3, 3, (mtx_type*)R);
}

void SixDOF::rotY(double phi, mtx_type* R)
{
    double c=cos(phi);
    double s=sin(phi);

    this->bufferR[0][0] = c;   this->bufferR[0][1] = 0; this->bufferR[0][2] = s;
    this->bufferR[1][0] = 0;   this->bufferR[1][1] = 1; this->bufferR[1][2] = 0;
    this->bufferR[2][0] = -s;  this->bufferR[2][1] = 0; this->bufferR[2][2] = c;
    
    Matrix.Multiply((mtx_type*)R, (mtx_type*)this->bufferR, 3, 3, 3, (mtx_type*)this->bufferR2);
    Matrix.Copy((mtx_type*)this->bufferR2, 3, 3, (mtx_type*)R);
}

void SixDOF::rotX(double psi, mtx_type* R)
{
    double c=cos(psi);
    double s=sin(psi);

    this->bufferR[0][0] = 1; this->bufferR[0][1] = 0; this->bufferR[0][2] = 0;
    this->bufferR[1][0] = 0; this->bufferR[1][1] = c; this->bufferR[1][2] = -s;
    this->bufferR[2][0] = 0; this->bufferR[2][1] = s; this->bufferR[2][2] = c;
    
    Matrix.Multiply((mtx_type*)R, (mtx_type*)this->bufferR, 3, 3, 3, (mtx_type*)this->bufferR2);
    Matrix.Copy((mtx_type*)this->bufferR2, 3, 3, (mtx_type*)R);
}


void SixDOF::rot2euler(mtx_type* R, mtx_type* euler) 
{
    // Method as per Paul, p 69.
    // euler = [phi theta psi]
    if (fabs(R[3 * 0 + 2]) < SINGULARITY_THRESHOLD && fabs(R[3 * 1 + 2]) < SINGULARITY_THRESHOLD)
    {
      // Singularity
      euler[0] = 0;
      double sp = 0;
      double cp = 1;
      euler[1] = atan2(cp*R[3 * 0 + 2] + sp*R[3 * 1 + 2], R[3 * 2 + 2]);
      euler[2] = atan2(-sp*R[3 * 0 + 0] + cp*R[3 * 1 + 0], -sp*R[3 * 0 + 1] + cp*R[3 * 1 + 1]);
    }
    else
    {
      // Non-singular
      euler[0] = atan2(R[3 * 1 + 2], R[3 * 2 + 2]);
      double sp = sin(this->pose[3]);
      double cp = cos(this->pose[3]);
      euler[1] = atan2(cp*R[3 * 0 + 2] + sp*R[3 * 1 + 2], R[3 * 2 + 2]);
      euler[2] = atan2(-sp*R[3 * 0 + 0] + cp*R[3 * 1 + 0], -sp*R[3 * 0 + 1] + cp*R[3 * 1 + 1]);
    }
}

void SixDOF::euler2rot(mtx_type* euler, mtx_type* R) 
{
    R[3 * 0 + 0] = 1; R[3 * 0 + 1] = 0; R[3 * 0 + 2] = 0;
    R[3 * 1 + 0] = 0; R[3 * 1 + 1] = 1; R[3 * 1 + 2] = 0;
    R[3 * 2 + 0] = 0; R[3 * 2 + 1] = 0; R[3 * 2 + 2] = 1;

    rotZ(euler[0], R);
    rotY(euler[1], R);
    rotZ(euler[2], R);
}

void SixDOF::euler2angvel(mtx_type* eulerFrom, mtx_type* eulerTo, mtx_type* angVel) 
{
  euler2rot((mtx_type*)eulerFrom,(mtx_type*)this->bufferR3);
  euler2rot((mtx_type*)eulerTo,(mtx_type*)this->bufferR4);
  Matrix.Copy((mtx_type*)this->bufferR3, 3, 3, (mtx_type*)this->bufferR2);

  Matrix.Transpose((mtx_type*)this->bufferR2, 3, 3, (mtx_type*)this->bufferR5);
  Matrix.Multiply((mtx_type*)this->bufferR4, (mtx_type*)this->bufferR5, 3, 3, 3, (mtx_type*)this->bufferR);
  rot2euler((mtx_type*)this->bufferR3, (mtx_type*)this->bufferV);
  rot2euler((mtx_type*)this->bufferR4, (mtx_type*)this->bufferU);
  rot2euler((mtx_type*)this->bufferR, (mtx_type*)this->bufferW);

  double c_z1 = cos(this->bufferV[0]);
  double s_z1 = sin(this->bufferV[0]);
  double c_y = sin(this->bufferV[1]);
  double s_y = sin(this->bufferV[1]);

  this->bufferR2[0][0] = 0; this->bufferR2[0][1] = -s_z1; this->bufferR2[0][2] = c_z1*s_y;
  this->bufferR2[1][0] = 0; this->bufferR2[1][1] = c_z1;  this->bufferR2[1][2] = s_z1*c_y;
  this->bufferR2[2][0] = 1; this->bufferR2[2][1] = 0;     this->bufferR2[2][2] = c_y;

  this->bufferW[0] = -((this->bufferV[0]+this->bufferW[0])-(this->bufferU[0]+this->bufferW[0]));
  this->bufferW[1] = -((this->bufferV[1]+this->bufferW[1])-(this->bufferU[1]+this->bufferW[1]));
  this->bufferW[2] = -((this->bufferV[0]+this->bufferW[2])-(this->bufferU[2]+this->bufferW[2]));

  Matrix.Multiply((mtx_type*)this->bufferR2, (mtx_type*)this->bufferW, 3, 3, 1, (mtx_type*)angVel);
}

//////////// ATTEMPT AT USING AXIS/ANGLE VECTORS INSTEAD OF EULER ANGLES /////////////////

/*
void SixDOF::rot2axis(mtx_type* R, mtx_type* ax)
{
    ax[0] = (R[3 * 2 + 1]-R[3 * 1 + 2]);
    ax[1] = (R[3 * 0 + 2]-R[3 * 2 + 0]);
    ax[3] = (R[3 * 1 + 0]-R[3 * 0 + 1]);
    ax[4] = 0;
    double norm = sqrt(ax[0]*ax[0] + ax[1]*ax[1] + ax[2]*ax[2]);
    double theta = atan2(0.5*norm,0.5*(R[3 * 0 + 0]+R[3 * 1 + 1]+R[3 * 2 + 2]-1));

    if (norm < SINGULARITY_THRESHOLD) {
      eigDecomp((mtx_type*)R, 3, 3, 10, (mtx_type*)this->bufferV);
      ax[0] = this->bufferV[0];
      ax[1] = this->bufferV[1];
      ax[2] = this->bufferV[2];
      norm = sqrt(ax[0]*ax[0] + ax[1]*ax[1] + ax[2]*ax[2]);
    }

    ax[0] /= norm;
    ax[1] /= norm;
    ax[3] /= norm;
    ax[4] = theta;
}

void SixDOF::axis2rot(mtx_type* ax, mtx_type* R)
{
    double c = cos(ax[3]);
    double s = sin(ax[3]);
    double t = 1 - c;
    double x = ax[0];
    double y = ax[1];
    double z = ax[2];
  
    R[3 * 0 + 0] = t*x*x + c;   R[3 * 0 + 1] = t*x*y - z*s; R[3 * 0 + 2] = t*x*z + y*s;
    R[3 * 1 + 0] = t*x*y + z*s; R[3 * 1 + 1] = t*y*y + c;   R[3 * 1 + 2] = t*y*z - x*s;
    R[3 * 2 + 0] = t*x*z - y*s; R[3 * 2 + 1] = t*y*z + x*s; R[3 * 2 + 2] = t*z*z + c;
}

// Power Itera
void SixDOF::eigDecomp(mtx_type* A, int m, int n, int iterations, mtx_type* out) 
{
  // Ideally choose a random vector
  // To decrease the chance that our vector
  // Is orthogonal to the eigenvector
  mtx_type* b_k = (mtx_type*)malloc(sizeof(mtx_type)*n);
  mtx_type* b_k1 = (mtx_type*)malloc(sizeof(mtx_type)*n);
  for (int i = 0; i < n; i++) {
    b_k[i] = ((double)random(100));
  }

  for (int i = 0; i < iterations; i++) {
    // calculate the matrix-by-vector product Ab
    Matrix.Multiply((mtx_type*)A, (mtx_type*)b_k, m, n, 1, (mtx_type*)b_k1);
    // calculate the norm
    double norm = 0;
    for (int j = 0; j < n; j++) {
      norm += (double)(b_k1[j]*b_k1[j]);
    }
    norm = sqrt(norm);
    // re-normalize the vector
    Matrix.Scale((mtx_type*)b_k1, n, 1, 1/norm);
    Matrix.Copy((mtx_type*)b_k1, n, 1, (mtx_type*)b_k);
    
  }
  
  Matrix.Copy((mtx_type*)b_k, n, 1, (mtx_type*)out);

  free(b_k);
  free(b_k1);
}
*/


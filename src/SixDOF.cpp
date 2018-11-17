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

#include "MemoryFree.h"

SixDOF::SixDOF(double* thetas, double* ds, double* alphas, double* as, int len)
{
    this->numJoints = len;
    this->numLinks = this->numJoints+1;
    
    this->links = new Link[this->numLinks];

    for (int i = 0; i < this->numLinks-1; i++)
    {
        this->links[i].fromJoint = new Joint[1];
        
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
            this->links[i-1].toJoint = new Joint[1];
            this->links[i-1].toJoint = links[i].fromJoint;
        }
    }

    this->links[this->numLinks-2].fromJoint = links[this->numLinks-3].toJoint;
    this->links[this->numLinks-2].isEndEffector = true;
    
    this->eye[0][0] = 1; this->eye[0][1] = 0; this->eye[0][2] = 0; this->eye[0][3] = 0;
    this->eye[1][0] = 0; this->eye[1][1] = 1; this->eye[1][2] = 0; this->eye[1][3] = 0;
    this->eye[2][0] = 0; this->eye[2][1] = 0; this->eye[2][2] = 1; this->eye[2][3] = 0;
    this->eye[3][0] = 0; this->eye[3][1] = 0; this->eye[3][2] = 0; this->eye[3][3] = 1;

    this->eye3[0][0] = 1; this->eye3[0][1] = 0; this->eye3[0][2] = 0;
    this->eye3[1][0] = 0; this->eye3[1][1] = 1; this->eye3[1][2] = 0;
    this->eye3[2][0] = 0; this->eye3[2][1] = 0; this->eye3[2][2] = 1;


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
    returnPose[0] = (double)this->pose[0];
    returnPose[1] = (double)this->pose[1];
    returnPose[2] = (double)this->pose[2];
    returnPose[3] = (double)this->pose[3];
    returnPose[4] = (double)this->pose[4];
    returnPose[5] = (double)this->pose[5];
}

void SixDOF::getJointAngles(double* returnAngles)
{
    for (int i = 0; i < this->numJoints; i++)
    {
        returnAngles[i] = this->links[i].fromJoint->a;
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

    Link endEffector = this->links[this->numJoints-1];
    
    this->pose[0] = endEffector.H0[0][3];
    this->pose[1] = endEffector.H0[1][3];
    this->pose[2] = endEffector.H0[2][3];

    mtx_type* bufferR = new mtx_type[3*3];
    mtx_type* bufferV = new mtx_type[3];

    bufferR[3*0+0] = endEffector.H0[0][0]; bufferR[3*0+1] = endEffector.H0[0][1]; bufferR[3*0+2] = endEffector.H0[0][2];
    bufferR[3*1+0] = endEffector.H0[1][0]; bufferR[3*1+1] = endEffector.H0[1][1]; bufferR[3*1+2] = endEffector.H0[1][2];
    bufferR[3*2+0] = endEffector.H0[2][0]; bufferR[3*2+1] = endEffector.H0[2][1]; bufferR[3*2+2] = endEffector.H0[2][2];
    
    rot2euler((mtx_type*)bufferR, (mtx_type*)bufferV);
    this->pose[3] = bufferV[0];
    this->pose[4] = bufferV[1];
    this->pose[5] = bufferV[2];

    delete[] bufferR;
    delete[] bufferV;
}



SixDOF::IKState_t SixDOF::inverseKinematics(double x, double y, double z, double phi, double theta, double psi, double grippingOffset, double alpha)
{

    bool failed       = false;
    double* dAngles   = new double[this->numJoints];
    mtx_type* bufferV = new mtx_type[3];
    mtx_type* bufferU = new mtx_type[3];
    mtx_type* bufferW = new mtx_type[3];
    mtx_type* bufferR = new mtx_type[3*3];
    mtx_type* oc      = new mtx_type[3];

    for (int i = 0; i < this->numJoints; i++){ dAngles[i] = 0; }
        
    int wristCenter = this->numJoints-1;
    for (int i = this->numJoints-2; i > 0; i--)
    {
        if (this->links[i].length == 0) {
          wristCenter = i; 
          break;
        }
    }

    int hasThirdOrientationJoint = 0;
    for (int i = wristCenter-1; i > 0; i--)
    {
        if (this->links[i].length == 0) {
          hasThirdOrientationJoint = 1; 
          break;
        }
    }


    bufferV[0] = x; 
    bufferV[1] = y; 
    bufferV[2] = z; 

    bufferU[0] = 0; 
    bufferU[1] = 0; 
    bufferU[2] = 1; 

    bufferW[0] = phi; 
    bufferW[1] = theta; 
    bufferW[2] = psi; 

    euler2rot((mtx_type*)bufferW, (mtx_type*)bufferR);
    //Matrix.Print((mtx_type*)bufferR, 3, 3, "bufferR");

    Matrix.Multiply((mtx_type*)bufferR, (mtx_type*)bufferU, 3, 3, 1, (mtx_type*)bufferW);
    Matrix.Scale((mtx_type*)bufferW, 3, 1, this->links[this->numJoints-1].length+grippingOffset);
    Matrix.Subtract((mtx_type*)bufferV, (mtx_type*)bufferW, 3, 1, (mtx_type*)oc);

    Matrix.Multiply((mtx_type*)bufferR, (mtx_type*)bufferU, 3, 3, 1, (mtx_type*)bufferW);
    Matrix.Scale((mtx_type*)bufferW, 3, 1, grippingOffset);
    Matrix.Subtract((mtx_type*)bufferV, (mtx_type*)bufferW, 3, 1, (mtx_type*)bufferU);

    x = bufferU[0];
    y = bufferU[1];
    z = bufferU[2];

    //Matrix.Print((mtx_type*)oc, 3, 1, "oc");

    int numSolve = wristCenter;
    
    // USE ITERATIVE METHOD TO SOLVE WRIST POSITION

    bufferU[0] = this->links[wristCenter-1].fromJoint->x;
    bufferU[1] = this->links[wristCenter-1].fromJoint->y;
    bufferU[2] = this->links[wristCenter-1].fromJoint->z;
    
    // Wrist center to Dest
    double b = vecDistance((mtx_type*)bufferU, (mtx_type*)oc);
    
    int* positionUpdateDir = new int[numSolve-1];
    for (int i = 0; i < numSolve; i++) {
        positionUpdateDir[i] = -b/fabs(b);
    }

    Serial.println("POSITION ---------------------");
    unsigned long time = millis();
    
    double prevB = b;
    int index = 0;
    int count = 0;
    double alphaTemp = alpha;
    while (b > REACHED_DESTINATION_THRESHOLD && count++ < INVKIN_TIMEOUT) 
    {
        //Serial.print(b);
        //Serial.print('\t');
        //Serial.println(alpha);
        //gamma = (a*a+b*b+c*c)/(2*a*b);
        
        double randVal = ((double)random(100))/10000;
        
        dAngles[index] = dAngles[index] + (alpha*b*b + randVal)*positionUpdateDir[index];
        alpha = alpha*0.999;

        forwardKinematics(dAngles);
        
        bufferU[0] = this->links[wristCenter-1].fromJoint->x;
        bufferU[1] = this->links[wristCenter-1].fromJoint->y;
        bufferU[2] = this->links[wristCenter-1].fromJoint->z;

        //Matrix.Print((mtx_type*)this->bufferU, 3, 1, "bufferU");
            
        // End effector to Dest
        b = vecDistance((mtx_type*)bufferU, (mtx_type*)oc);

        if (b <= prevB) {
            if (index == numSolve-hasThirdOrientationJoint-1) 
                index = 0;
            else
                index++;
        } else {
            positionUpdateDir[index] = -positionUpdateDir[index];
        }
        
        prevB = b;
    }

    if (count > INVKIN_TIMEOUT) failed = true;
    
    forwardKinematics((double*)dAngles);

    delete[] positionUpdateDir;

    time = millis() - time;
    Serial.print("TIME ELAPSED: ");Serial.print(time);Serial.println(" ms");
    Serial.print("UNCERTAINTY: ");Serial.println(b);
    
    ///////////////////////
    // USE ITERATIVE METHOD TO SOLVE ORIENTATION

    bufferU[0] = x; 
    bufferU[1] = y; 
    bufferU[2] = z; 
    
    bufferV[0] = this->pose[0];
    bufferV[1] = this->pose[1];
    bufferV[2] = this->pose[2];
    
    // End effector to Dest
    b = vecDistance((mtx_type*)bufferV, (mtx_type*)bufferU);
    
    int* orientationUpdateDir = new int[this->numJoints-numSolve-1];
    for (int i = 0; i < this->numJoints-numSolve; i++) {
        orientationUpdateDir[i] = b/fabs(b);
    }

    Serial.println("ORIENTATION ---------------------");
    time = millis();

    prevB = b;
    index = numSolve-hasThirdOrientationJoint;
    count = 0;
    alpha = alphaTemp;
    while (b > REACHED_DESTINATION_THRESHOLD && count++ < INVKIN_TIMEOUT) 
    {
        //Serial.print(b);
        //Serial.print('\t');
        //Serial.println(alpha);
        //gamma = (a*a+b*b+c*c)/(2*a*b);
        dAngles[index] = dAngles[index] + (alpha*fabs(b))*orientationUpdateDir[index-numSolve];
        alpha = alpha*0.999;


        forwardKinematics((double*)dAngles);
        
        bufferV[0] = this->pose[0];
        bufferV[1] = this->pose[1];
        bufferV[2] = this->pose[2];
            
        // End effector to Dest
        b = vecDistance((mtx_type*)bufferV, (mtx_type*)bufferU);
    
        if (b <= prevB) {
            if (index == this->numJoints-1) 
                index = numSolve-hasThirdOrientationJoint;
            else
                index++;
        } else {
            orientationUpdateDir[index-numSolve] = -orientationUpdateDir[index-numSolve];
        }
        
        prevB = b;
    }

    if (count > INVKIN_TIMEOUT) failed = true;
    
    forwardKinematics((double*)dAngles);
    time = millis() - time;
    Serial.print("TIME ELAPSED: ");Serial.print(time);Serial.println(" ms");
    Serial.print("UNCERTAINTY: ");Serial.println(b);

    Serial.println("ORIENTATION ADJUST ---------------------");
    time = millis();

    int timeout = 0;
    double optimalAngle = 0;
    double bestError = 1000;
    do {
      bufferU[0] = this->links[this->numJoints-1].H0[0][1];
      bufferU[1] = this->links[this->numJoints-1].H0[1][1];
      bufferU[2] = this->links[this->numJoints-1].H0[2][1];
  
      bufferV[0] = bufferR[3*0+1]; 
      bufferV[1] = bufferR[3*1+1];
      bufferV[2] = bufferR[3*2+1]; 
  
      double sU = norm(bufferU);
      double sV = norm(bufferV);
      double cU = dot(bufferU, bufferV);
      double error = acos(cU/(sV*sU));
  
      dAngles[this->numJoints-1] = dAngles[this->numJoints-1]-error*0.5;
      if (error < bestError) {
        bestError = error;
        optimalAngle = dAngles[this->numJoints-1];
        timeout = 0;
      } else {
        timeout++;
      }
      forwardKinematics((double*)dAngles);
      //Serial.println(bestError);
    } while (timeout < 10);

    dAngles[this->numJoints-1] = optimalAngle;
    
    for (int i = 0; i < this->numJoints; i++) dAngles[i] = fmod(dAngles[i], 2*PI);
    
    forwardKinematics((double*)dAngles);
    time = millis() - time;
    Serial.print("TIME ELAPSED: ");Serial.print(time);Serial.println(" ms");
    Serial.print("UNCERTAINTY: ");Serial.println(bestError);

    delete[] orientationUpdateDir;
    delete[] dAngles;
    delete[] bufferV;
    delete[] bufferU;
    delete[] bufferW;
    delete[] bufferR;
    delete[] oc;

    if (failed == false) {
      Serial.println("SUCCESS!");
      this->ikStatus = SUCCESS;
      return SUCCESS;
    } else {
      Serial.println("FAIL!");
      this->ikStatus = FAILED;
      return FAILED;
    }
}

SixDOF::IKState_t SixDOF::inverseKinematics2(double x, double y, double z, double phi, double theta, double psi, double grippingOffset, double alpha)
{
    mtx_type* dAngles   = new mtx_type[this->numJoints];
    mtx_type* error = new mtx_type[6];
    mtx_type* bufferAxis = new mtx_type[4];
    mtx_type* bufferV = new mtx_type[3];
    mtx_type* bufferR = new mtx_type[3*3];
    mtx_type* bufferR2 = new mtx_type[3*3];
    mtx_type* bufferR3 = new mtx_type[3*3];
    mtx_type* oc = new mtx_type[3];
    mtx_type* J = new mtx_type[this->numJoints*6];
    mtx_type* transposeJ = new mtx_type[this->numJoints*6];
    mtx_type* bufferinvJ = new mtx_type[this->numJoints*this->numJoints];
    mtx_type* invJ = new mtx_type[this->numJoints*6];

    for (int i = 0; i < this->numJoints; i++) dAngles[i] = 0;

    bufferV[0] = phi; 
    bufferV[1] = theta; 
    bufferV[2] = psi; 

    euler2rot((mtx_type*)bufferV, (mtx_type*)bufferR);

    this->ikStatus = SUCCESS;

    while(true)
    {

      Link endEffector = this->links[this->numJoints-1];

      error[0] = x-endEffector.H0[0][3];
      error[1] = y-endEffector.H0[1][3];
      error[2] = z-endEffector.H0[2][3];

      t2r((mtx_type*)this->links[this->numJoints-1].H0,(mtx_type*)bufferR2);

      Matrix.Transpose((mtx_type*)bufferR2, 3, 3, (mtx_type*)bufferR3);
      Matrix.Multiply((mtx_type*)bufferR3, (mtx_type*)bufferR, 3, 3, 3, (mtx_type*)bufferR2);
      
      rot2axis(bufferR2, bufferAxis);
      error[3] = bufferAxis[0]*bufferAxis[3];
      error[4] = bufferAxis[1]*bufferAxis[3];
      error[5] = bufferAxis[2]*bufferAxis[3];

      computeJointJacobian((mtx_type*)endEffector.H0, (mtx_type*)this->eye, (mtx_type*)this->links[0].fromJoint->J);
      J[this->numJoints * 0 + 0] = this->links[0].fromJoint->J[0];
      J[this->numJoints * 1 + 0] = this->links[0].fromJoint->J[1];
      J[this->numJoints * 2 + 0] = this->links[0].fromJoint->J[2];
      J[this->numJoints * 3 + 0] = this->links[0].fromJoint->J[3];
      J[this->numJoints * 4 + 0] = this->links[0].fromJoint->J[4];
      J[this->numJoints * 5 + 0] = this->links[0].fromJoint->J[5];
      
      for (int i = 1; i < this->numJoints; i++)
      {
          computeJointJacobian((mtx_type*)endEffector.H0, (mtx_type*)this->links[i-1].H0, (mtx_type*)this->links[i].fromJoint->J);
          
          J[this->numJoints * 0 + i] = this->links[i].fromJoint->J[0];
          J[this->numJoints * 1 + i] = this->links[i].fromJoint->J[1];
          J[this->numJoints * 2 + i] = this->links[i].fromJoint->J[2];
          J[this->numJoints * 3 + i] = this->links[i].fromJoint->J[3];
          J[this->numJoints * 4 + i] = this->links[i].fromJoint->J[4];
          J[this->numJoints * 5 + i] = this->links[i].fromJoint->J[5];
      }

      //Matrix.Print((mtx_type*)J, 6, this->numJoints, "J");

      Matrix.Transpose((mtx_type*)J, 6, this->numJoints, (mtx_type*)transposeJ);
      Matrix.Multiply((mtx_type*)transposeJ, (mtx_type*)J, this->numJoints, 6, this->numJoints, (mtx_type*)bufferinvJ);
      
      int success = Matrix.Invert((mtx_type*)bufferinvJ, this->numJoints);
      if (success == 0) {
        this->ikStatus = FAILED;
        break;
      }

      Matrix.Multiply((mtx_type*)bufferinvJ, (mtx_type*)transposeJ, this->numJoints, this->numJoints, 6, (mtx_type*)invJ);

      Matrix.Multiply((mtx_type*)invJ, (mtx_type*)error, this->numJoints, 6, 1, (mtx_type*)dAngles);

      for (int i = 0; i < this->numJoints; i++)
      {
          dAngles[i] = this->links[i].fromJoint->a + alpha*dAngles[i];
          dAngles[i] = fmod(dAngles[i], 2*PI);
      }

      alpha = max(0.01, alpha*0.99);

      forwardKinematics((double*)dAngles);

      double np = sqrt(error[0]*error[0] + error[1]*error[1] + error[2]*error[2]);
      double no = sqrt(error[3]*error[3] + error[4]*error[4] + error[5]*error[5]);
      Serial.print(alpha);
      Serial.print("\t");
      Serial.print(np);
      Serial.print("\t");
      Serial.println(no);

      if (np < REACHED_DESTINATION_THRESHOLD && no < REACHED_DESTINATION_THRESHOLD) break;

    }

    delete[] error;
    delete[] dAngles;
    delete[] bufferAxis;
    delete[] bufferV;
    delete[] bufferR;
    delete[] bufferR2;
    delete[] bufferR3;
    delete[] oc;
    delete[] J;
    delete[] transposeJ;
    delete[] bufferinvJ;
    delete[] invJ;

    if (this->ikStatus == SUCCESS) {
      Serial.println("SUCCESS!");
      this->ikStatus = SUCCESS;
      return SUCCESS;
    } else {
      Serial.println("FAIL!");
      this->ikStatus = FAILED;
      return FAILED;
    }
}

//////////////////////////////////////////

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
    mtx_type* zPrev = new mtx_type[3];
    mtx_type* oPrev = new mtx_type[3];
    mtx_type* oEnd  = new mtx_type[3];
    mtx_type* bufferV  = new mtx_type[3];
    mtx_type* bufferU  = new mtx_type[3];
    
    zPrev[0] = prevH[4 * 0 + 2];
    zPrev[1] = prevH[4 * 1 + 2];
    zPrev[2] = prevH[4 * 2 + 2];
    
    oPrev[0] = prevH[4 * 0 + 3];
    oPrev[1] = prevH[4 * 1 + 3];
    oPrev[2] = prevH[4 * 2 + 3];
    
    oEnd[0] = endH[4 * 0 + 3];
    oEnd[1] = endH[4 * 1 + 3];
    oEnd[2] = endH[4 * 2 + 3];

    Matrix.Subtract((mtx_type*)oEnd, (mtx_type*)oPrev, 3, 1, (mtx_type*)bufferU);
    crossProduct(zPrev, bufferU, bufferV);
    J[0] = bufferV[0];
    J[1] = bufferV[1];
    J[2] = bufferV[2];
    J[3] = zPrev[0];
    J[4] = zPrev[1];
    J[5] = zPrev[2];

    delete[] bufferU;
    delete[] bufferV;
    delete[] zPrev;
    delete[] oPrev;
    delete[] oEnd;
}

void SixDOF::screwZ(double theta, double d, mtx_type* H)
{
    double c=cos(theta);
    double s=sin(theta);

    mtx_type* bufferH = new mtx_type[4*4];
    mtx_type* bufferH2 = new mtx_type[4*4];

    bufferH[4*0+0] = c;  bufferH[4*0+1] = -s; bufferH[4*0+2] = 0; bufferH[4*0+3] = 0;
    bufferH[4*1+0] = s;  bufferH[4*1+1] = c;  bufferH[4*1+2] = 0; bufferH[4*1+3] = 0;
    bufferH[4*2+0] = 0;  bufferH[4*2+1] = 0;  bufferH[4*2+2] = 1; bufferH[4*2+3] = d;
    bufferH[4*3+0] = 0;  bufferH[4*3+1] = 0;  bufferH[4*3+2] = 0; bufferH[4*3+3] = 1;
    
    Matrix.Multiply((mtx_type*)H, (mtx_type*)bufferH, 4, 4, 4, (mtx_type*)bufferH2);
    Matrix.Copy((mtx_type*)bufferH2, 4, 4, (mtx_type*)H);
    delete[] bufferH;
    delete[] bufferH2;
}

void SixDOF::screwX(double alpha, double a, mtx_type* H)
{
    double c=cos(alpha);
    double s=sin(alpha);

    mtx_type* bufferH = new mtx_type[4*4];
    mtx_type* bufferH2 = new mtx_type[4*4];

    bufferH[4*0+0] = 1;  bufferH[4*0+1] = 0;  bufferH[4*0+2] = 0;  bufferH[4*0+3] = a;
    bufferH[4*1+0] = 0;  bufferH[4*1+1] = c;  bufferH[4*1+2] = -s; bufferH[4*1+3] = 0;
    bufferH[4*2+0] = 0;  bufferH[4*2+1] = s;  bufferH[4*2+2] = c;  bufferH[4*2+3] = 0;
    bufferH[4*3+0] = 0;  bufferH[4*3+1] = 0;  bufferH[4*3+2] = 0;  bufferH[4*3+3] = 1;
    
    Matrix.Multiply((mtx_type*)H, (mtx_type*)bufferH, 4, 4, 4, (mtx_type*)bufferH2);
    Matrix.Copy((mtx_type*)bufferH2, 4, 4, (mtx_type*)H);
    delete[] bufferH;
    delete[] bufferH2;
}

void SixDOF::crossProduct(mtx_type* v, mtx_type* u, mtx_type* w) 
{
    double jk = v[1]*u[2]-u[1]*v[2];
    double ik = v[0]*u[2]-u[0]*v[2];
    double ij = v[0]*u[1]-u[0]*v[1];
    
    Matrix.Scale((mtx_type*)this->v_i, 3, 1, jk);
    Matrix.Scale((mtx_type*)this->v_j, 3, 1, ik);
    Matrix.Scale((mtx_type*)this->v_k, 3, 1, ij);

    mtx_type* bufferV = new mtx_type[3];
    
    Matrix.Subtract((mtx_type*)this->v_i, (mtx_type*)this->v_j, 3, 1, (mtx_type*)bufferV);
    Matrix.Add((mtx_type*)bufferV, (mtx_type*)this->v_k, 3, 1, (mtx_type*)w);

    this->v_i[0] = 1;
    this->v_i[1] = 0;
    this->v_i[2] = 0;
    this->v_j[0] = 0;
    this->v_j[1] = 1;
    this->v_j[2] = 0;
    this->v_k[0] = 0;
    this->v_k[1] = 0;
    this->v_k[2] = 1;

    delete[] bufferV;
}


double SixDOF::angleDiff(mtx_type* v, mtx_type* u)
{
    /*double dotProd = v[0]*u[0] + v[1]*u[1] + v[2]*u[2];
    double normV = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    double normU = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    
    double c_ang = dotProd/(normV*normU);
    double s_ang = sqrt(1-c_ang*c_ang);
    
    return atan2(s_ang, c_ang);*/

    mtx_type* w = new mtx_type[3];
    crossProduct((mtx_type*) v, (mtx_type*) u, (mtx_type*) w);
    double normW = sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
    double normU = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    double normV = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    delete[] w;
    return asin(normW/(normU*normV));
} 


double SixDOF::vecDistance(mtx_type* from, mtx_type* to)
{
  double diffX = to[0]-from[0];
  double diffY = to[1]-from[1];
  double diffZ = to[2]-from[2];
  return sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ);
}

double SixDOF::norm(mtx_type* v)
{
  double x = v[0];
  double y = v[1];
  double z = v[2];
  return sqrt(x*x + y*y + z*z);
}
/*
double SixDOF::norm6(mtx_type* v)
{
  double x = v[0];
  double y = v[1];
  double z = v[2];
  double psi = v[3];
  double phi = v[4];
  double the = v[5];
  return sqrt(x*x + y*y + z*z + psi*psi + phi*phi + the*the);
}*/

double SixDOF::dot(mtx_type* u, mtx_type* v)
{
  double x = u[0]*v[0];
  double y = u[1]*v[1];
  double z = u[2]*v[2];
  return x+y+z;
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
      mtx_type* temp_M = new mtx_type[splitSize * splitSize];
      
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
      delete[] temp_M;
  }
  return det;
}

void SixDOF::rotZ(double theta, mtx_type* R)
{
    double c=cos(theta);
    double s=sin(theta);

    mtx_type* bufferR = new mtx_type[3*3];
    mtx_type* result  = new mtx_type[3*3];

    bufferR[3*0+0] = c;  bufferR[3*0+1] = -s; bufferR[3*0+2] = 0;
    bufferR[3*1+0] = s;  bufferR[3*1+1] = c;  bufferR[3*1+2] = 0;
    bufferR[3*2+0] = 0;  bufferR[3*2+1] = 0;  bufferR[3*2+2] = 1;
    
    Matrix.Multiply((mtx_type*)R, (mtx_type*)bufferR, 3, 3, 3, (mtx_type*)result);
    Matrix.Copy((mtx_type*)result, 3, 3, (mtx_type*)R);
    delete[] bufferR;
    delete[] result;
}

void SixDOF::rotY(double phi, mtx_type* R)
{
    double c=cos(phi);
    double s=sin(phi);

    mtx_type* bufferR = new mtx_type[3*3];
    mtx_type* result  = new mtx_type[3*3];

    bufferR[3*0+0] = c;  bufferR[3*0+1] = 0; bufferR[3*0+2] = s;
    bufferR[3*1+0] = 0;  bufferR[3*1+1] = 1; bufferR[3*1+2] = 0;
    bufferR[3*2+0] = -s; bufferR[3*2+1] = 0; bufferR[3*2+2] = c;
    
    Matrix.Multiply((mtx_type*)R, (mtx_type*)bufferR, 3, 3, 3, (mtx_type*)result);
    Matrix.Copy((mtx_type*)result, 3, 3, (mtx_type*)R);
    delete[] bufferR;
    delete[] result;
}

void SixDOF::rotX(double psi, mtx_type* R)
{
    double c=cos(psi);
    double s=sin(psi);

    mtx_type* bufferR = new mtx_type[3*3];
    mtx_type* result  = new mtx_type[3*3];

    bufferR[3*0+0] = 1; bufferR[3*0+1] = 0; bufferR[3*0+2] = 0;
    bufferR[3*1+0] = 0; bufferR[3*1+1] = c; bufferR[3*1+2] = -s;
    bufferR[3*2+0] = 0; bufferR[3*2+1] = s; bufferR[3*2+2] = c;
    
    Matrix.Multiply((mtx_type*)R, (mtx_type*)bufferR, 3, 3, 3, (mtx_type*)result);
    Matrix.Copy((mtx_type*)result, 3, 3, (mtx_type*)R);
    delete[] bufferR;
    delete[] result;
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


void SixDOF::t2r(mtx_type* H, mtx_type* R) 
{
    R[3 * 0 + 0] = H[4 * 0 + 0]; R[3 * 0 + 1] = H[4 * 0 + 1]; R[3 * 0 + 2] = H[4 * 0 + 2];
    R[3 * 1 + 0] = H[4 * 1 + 0]; R[3 * 1 + 1] = H[4 * 1 + 1]; R[3 * 1 + 2] = H[4 * 1 + 2];
    R[3 * 2 + 0] = H[4 * 2 + 0]; R[3 * 2 + 1] = H[4 * 2 + 1]; R[3 * 2 + 2] = H[4 * 2 + 2];
}

double SixDOF::traceR(mtx_type* R) 
{
    return R[3 * 0 + 0]+R[3 * 1 + 1]+R[3 * 2 + 2];
}

double SixDOF::maxDiagR(mtx_type* R) 
{
    return max(R[3 * 0 + 0], max(R[3 * 1 + 1], R[3 * 2 + 2]));
}

int SixDOF::maxDiagRCol(mtx_type* R) 
{
    if (R[3 * 0 + 0] > R[3 * 1 + 1] && R[3 * 0 + 0] > R[3 * 2 + 2])
      return 0;
    if (R[3 * 1 + 1] > R[3 * 0 + 0] && R[3 * 1 + 1] > R[3 * 2 + 2])
      return 1;
    if (R[3 * 2 + 2] > R[3 * 1 + 1] && R[3 * 2 + 2] > R[3 * 0 + 0])
      return 2;
}

/*
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
}*/

//////////// ATTEMPT AT USING AXIS/ANGLE VECTORS INSTEAD OF EULER ANGLES /////////////////


void SixDOF::rot2axis(mtx_type* R, mtx_type* ax)
{
    if (fabs(traceR(R)-3) < SINGULARITY_THRESHOLD) {
      ax[0] = 0;
      ax[1] = 0;
      ax[2] = 0;
      ax[3] = 0;
    } else if (fabs(traceR(R)+1) < SINGULARITY_THRESHOLD) {
      double mx = maxDiagR(R);
      int mxCol = maxDiagRCol(R);

      double lg = sqrt(2*(1+mx));
      ax[0] = (R[3*0+mxCol] + this->eye3[0][mxCol])/lg;
      ax[1] = (R[3*1+mxCol] + this->eye3[1][mxCol])/lg;
      ax[2] = (R[3*2+mxCol] + this->eye3[2][mxCol])/lg;
      
      ax[3] = PI;
    } else {
      ax[3] = acos((traceR(R)-1)/2);

      mtx_type* skw = new mtx_type[3*3];
      mtx_type* transposeR = new mtx_type[3*3];

      Matrix.Transpose((mtx_type*)R, 3, 3, (mtx_type*)transposeR);
      Matrix.Subtract((mtx_type*)R, (mtx_type*)transposeR, 3, 3, (mtx_type*)skw);
      Matrix.Scale((mtx_type*)skw, 3, 3, 0.5/sin(ax[3]));
      
      ax[0] = 0.5*(skw[3*2+1]-skw[3*1+2]);
      ax[1] = 0.5*(skw[3*0+2]-skw[3*2+0]);
      ax[2] = 0.5*(skw[3*1+0]-skw[3*0+1]);

      delete[] skw;
      delete[] transposeR;
    }
}



/*
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
*/

// Power Itera
/*
void SixDOF::eigDecomp(mtx_type* A, int m, int n, int iterations, mtx_type* out) 
{
  // Ideally choose a random vector
  // To decrease the chance that our vector
  // Is orthogonal to the eigenvector
  mtx_type* b_k = new mtx_type[n];
  mtx_type* b_k1 = new mtx_type[n];
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

  delete[] b_k;
  delete[] b_k1;
}
*/


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


SixDOF::SixDOF(float* thetas, float* ds, float* alphas, float* as, int len)
{

    this->numJoints = len;
    this->numLinks = this->numJoints+1;
    
    this->links = (Link*)malloc(sizeof(Link)*(this->numLinks));

    for (int i = 0; i < this->numLinks; i++)
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
        
        this->links[i].length = ds[i];

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
}

void SixDOF::getPose(float* returnPose)
{
    returnPose[0] = (float)this->pose[0];
    returnPose[1] = (float)this->pose[1];
    returnPose[2] = (float)this->pose[2];
    returnPose[3] = (float)this->pose[3];
    returnPose[4] = (float)this->pose[4];
    returnPose[5] = (float)this->pose[5];
}

void SixDOF::getJointAngles(float* returnAngles)
{
    for (int i = 1; i < this->numLinks; i++)
    {
        returnAngles[i-1] = this->links[i-1].fromJoint->a;
    }
}

void SixDOF::forwardKinematics(float* angles)
{
    Matrix.Copy((mtx_type*)this->eye, 4, 4, (mtx_type*)this->links[0].H0);
    Matrix.Copy((mtx_type*)this->eye, 4, 4, (mtx_type*)this->links[0].fromJoint->H);
    float theta = this->links[0].fromJoint->DH_Params[0];
    float d     = this->links[0].fromJoint->DH_Params[1];
    float alpha = this->links[0].fromJoint->DH_Params[2];
    float a     = this->links[0].fromJoint->DH_Params[3];
    computeDH(angles[0]+theta, d, alpha, a, (mtx_type*)this->links[0].H0);
    computeDH(angles[0]+theta, d, alpha, a, (mtx_type*)this->links[0].fromJoint->H);
    //Matrix.Print((mtx_type*)this->links[0].H0, 4, 4, "H0");

    this->links[0].fromJoint->x = this->links[0].H0[0][3];
    this->links[0].fromJoint->y = this->links[0].H0[1][3];
    this->links[0].fromJoint->z = this->links[0].H0[2][3];
    
    this->links[0].fromJoint->a = angles[0];
    
    for (int i = 1; i < this->numLinks-1; i++)
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

        //Matrix.Print((mtx_type*)this->links[i].H0, 4, 4, "H");
    }
    
    Link endEffector = this->links[this->numLinks-2];
    
    this->pose[0] = endEffector.H0[0][3];
    this->pose[1] = endEffector.H0[1][3];
    this->pose[2] = endEffector.H0[2][3];
    
    this->bufferV[0] = this->links[0].H0[0][0];
    this->bufferV[1] = this->links[0].H0[1][0];
    this->bufferV[2] = this->links[0].H0[2][0];
    this->bufferU[0] = endEffector.H0[0][0];
    this->bufferU[1] = endEffector.H0[1][0];
    this->bufferU[2] = endEffector.H0[2][0];
    
    this->pose[3] = angleDiff(bufferV, bufferU);
    
    this->bufferV[0] = this->links[0].H0[0][1];
    this->bufferV[1] = this->links[0].H0[1][1];
    this->bufferV[2] = this->links[0].H0[2][1];
    this->bufferU[0] = endEffector.H0[0][1];
    this->bufferU[1] = endEffector.H0[1][1];
    this->bufferU[2] = endEffector.H0[2][1];
    
    this->pose[4] = angleDiff(this->bufferV, this->bufferU)+PI/2;
    
    this->bufferV[0] = this->links[0].H0[0][2];
    this->bufferV[1] = this->links[0].H0[1][2];
    this->bufferV[2] = this->links[0].H0[2][2];
    this->bufferU[0] = endEffector.H0[0][2];
    this->bufferU[1] = endEffector.H0[1][2];
    this->bufferU[2] = endEffector.H0[2][2];
    
    this->pose[5] = angleDiff(this->bufferV, this->bufferU);
}

int SixDOF::inverseKinematics(float x, float y, float z, float theta, float phi, float psi, float dt)
{
    Link endEffector = this->links[this->numLinks-2];
    
    this->newPose[0] = x;
    this->newPose[1] = y;
    this->newPose[2] = z;
    this->newPose[3] = theta;
    this->newPose[4] = phi;
    this->newPose[5] = psi;

    for(int i = 0; i < 6; i++)
          for(int j = 0; j < this->numJoints; j++)
                this->J[this->numJoints*i + j] = 0;

    ////Matrix.Print((mtx_type*)this->J, 6, this->numJoints, "J");
    
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

    //for(int i = 0; i < 6; i++)
    //      for(int j = 0; j < this->numJoints; j++)
    //            if (fabs(this->J[this->numJoints*i + j]) < 0.000001) this->J[this->numJoints*i + j] = 0;
    

    Matrix.Print((mtx_type*)this->J, 6, this->numJoints, "J");
    //Matrix.Print((mtx_type*)this->bufferinvJ, this->numJoints, this->numJoints, "J'J");
    
    // Compute determinant of J > JACOBIAN_SINGULARITY_THRESHOLD
    float detJ = determinant((mtx_type*)this->J, this->numJoints);
    if (abs(detJ) > JACOBIAN_SINGULARITY_THRESHOLD)
    {
        int success = 0;
#ifdef INVERSE_JACOBIAN_METHOD
        if (this->numJoints == 6) {
          Matrix.Copy((mtx_type*)this->J, this->numJoints, this->numJoints, (mtx_type*)this->invJ);
          success = Matrix.Invert((mtx_type*)this->invJ, this->numJoints);
          if (success == 0) return 0;
          //Matrix.Scale((mtx_type*)this->invJ, this->numJoints, this->numJoints, -1);
        } else {
          Matrix.Transpose((mtx_type*)this->J, 6, this->numJoints, (mtx_type*)this->transposeJ);
          Matrix.Multiply((mtx_type*)this->transposeJ, (mtx_type*)this->J, this->numJoints, 6, this->numJoints, (mtx_type*)this->bufferinvJ);
          
          success = Matrix.Invert((mtx_type*)this->bufferinvJ, this->numJoints);
          if (success == 0) return 0;
          //Matrix.Scale((mtx_type*)this->bufferinvJ, this->numJoints, this->numJoints, -1);

          Matrix.Multiply((mtx_type*)this->bufferinvJ, (mtx_type*)this->transposeJ, this->numJoints, this->numJoints, 6, (mtx_type*)this->invJ);
        }
#endif

#ifdef PSEUDO_INVERSE_JACOBIAN_METHOD
        Matrix.Transpose((mtx_type*)this->J, 6, this->numJoints, (mtx_type*)this->transposeJ);
        Matrix.Multiply((mtx_type*)this->transposeJ, (mtx_type*)this->J, this->numJoints, 6, this->numJoints, (mtx_type*)this->bufferinvJ);
        
        success = Matrix.Invert((mtx_type*)this->bufferinvJ, this->numJoints);
        if (success == 0) return 0;
        //Matrix.Scale((mtx_type*)this->bufferinvJ, this->numJoints, this->numJoints, -1);

        Matrix.Multiply((mtx_type*)this->bufferinvJ, (mtx_type*)this->transposeJ, this->numJoints, this->numJoints, 6, (mtx_type*)this->invJ);
#endif

        //Matrix.Print((mtx_type*)this->invJ, this->numJoints, 6, "invJ");

        //Matrix.Print((mtx_type*)this->newPose, 6, 1, "newPose");
        //Matrix.Print((mtx_type*)this->pose, 6, 1, "pose");
        
        Matrix.Subtract((mtx_type*)this->newPose, (mtx_type*)this->pose, 6, 1, (mtx_type*)this->bufferPose);
        //Matrix.Print((mtx_type*)this->bufferPose, 6, 1, "diffPose");
        Matrix.Copy((mtx_type*)this->bufferPose, 6, 1, (mtx_type*)this->newPose);
        Matrix.Scale((mtx_type*)this->newPose, 6, 1, ALPHA);
        Matrix.Multiply((mtx_type*)this->invJ, (mtx_type*)this->newPose, this->numJoints, 6, 1, (mtx_type*)this->dAngles);
        
        for (int i = 0; i < this->numLinks-1; i++)
        {
            //Serial.print(i);
            //Serial.print("\tOld: "); Serial.print(this->links[i].fromJoint->a*180.0f/PI);
            //if (abs(this->dAngles[i]) > 1000) this->dAngles[i] = 0;
            Serial.print("\tDiff: "); Serial.print(this->dAngles[i]*dt*180.0f/PI);
            this->dAngles[i] = this->links[i].fromJoint->a + this->dAngles[i];
            //Serial.print("\tNew: "); Serial.println(this->dAngles[i]*180.0f/PI);
            Serial.println();
        }
        
        forwardKinematics((float*)this->dAngles);
    } else { 
      Serial.println("Jacobian singular - Determinant below defined JACOBIAN_SINGULARITY_THRESHOLD");
      return 0; 
    }
    
    return 1;
}

void SixDOF::computeDH(float theta, float d, float alpha, float a, mtx_type* H) 
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

void SixDOF::screwZ(float theta, float d, mtx_type* H)
{
    float c=cos(theta);
    float s=sin(theta);

    this->bufferH[0][0] = c;  this->bufferH[0][1] = -s; this->bufferH[0][2] = 0; this->bufferH[0][3] = 0;
    this->bufferH[1][0] = s; this->bufferH[1][1] = c; this->bufferH[1][2] = 0; this->bufferH[1][3] = 0;
    this->bufferH[2][0] = 0;  this->bufferH[2][1] = 0; this->bufferH[2][2] = 1; this->bufferH[2][3] = d;
    this->bufferH[3][0] = 0;  this->bufferH[3][1] = 0; this->bufferH[3][2] = 0; this->bufferH[3][3] = 1;
    
    Matrix.Multiply((mtx_type*)H, (mtx_type*)this->bufferH, 4, 4, 4, (mtx_type*)this->bufferH2);
    Matrix.Copy((mtx_type*)this->bufferH2, 4, 4, (mtx_type*)H);
}

void SixDOF::screwX(float alpha, float a, mtx_type* H)
{
    float c=cos(alpha);
    float s=sin(alpha);

    this->bufferH[0][0] = 1; this->bufferH[0][1] = 0; this->bufferH[0][2] = 0;  this->bufferH[0][3] = a;
    this->bufferH[1][0] = 0; this->bufferH[1][1] = c; this->bufferH[1][2] = -s; this->bufferH[1][3] = 0;
    this->bufferH[2][0] = 0; this->bufferH[2][1] = s; this->bufferH[2][2] = c;  this->bufferH[2][3] = 0;
    this->bufferH[3][0] = 0; this->bufferH[3][1] = 0; this->bufferH[3][2] = 0;  this->bufferH[3][3] = 1;
    
    Matrix.Multiply((mtx_type*)H, (mtx_type*)this->bufferH, 4, 4, 4, (mtx_type*)this->bufferH2);
    Matrix.Copy((mtx_type*)this->bufferH2, 4, 4, (mtx_type*)H);
}

void SixDOF::crossProduct(mtx_type* v, mtx_type* u, mtx_type* w) 
{
    float jk = v[1]*u[2]-u[1]*v[2];
    float ik = v[0]*u[2]-u[0]*v[2];
    float ij = v[0]*u[1]-u[0]*v[1];
    
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

float SixDOF::angleDiff(mtx_type* v, mtx_type* u)
{
    float dotProd = v[0]*u[0] + v[1]*u[1] + v[2]*u[2];
    float normV = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    float normU = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    
    float c_ang = dotProd/(normV*normU);
    float s_ang = sqrt(1-c_ang*c_ang);
    
    return atan2(s_ang, c_ang);
} 

float SixDOF::determinant(mtx_type* M, int len) 
{
  float det = 0;
  if (len == 2) 
  {
      det = (float)(M[len * 0 + 0]*M[len * 1 + 1]-M[len * 0 + 1]*M[len * 1 + 0]);
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
              det = (float)(M[len * 0 + i]*o*determinant((mtx_type*)temp_M, splitSize));
          } else {
              det = det + (float)(M[len * 0 + i]*o*determinant((mtx_type*)temp_M, splitSize));
          }
      }
      free((mtx_type*)temp_M);
  }
  return det;
}

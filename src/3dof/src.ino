
#include "SixDOF.h"
#include "MatrixMath.h"
#include "MemoryFree.h"

double thetas[] = {0, PI/2, -PI/2, 0};
double ds[] = {0.25, 0, 0, 0.2};
double alphas[] = {PI/2, 0, -PI/2, 0};
double as[] = {0, 0.25, 0, 0};

double jointAngles[] = {0, PI/3, PI/2, 0};
double jointAngles2[] = {0, PI/4, 0, 0};

double pose[] = {0, 0, 0, 0, 0, 0};

double desiredPose1[] = {0, -0.25, 0.05, 0, 3.1416, 1.5708};

SixDOF manipulator(thetas, ds, alphas, as, 4);

void setup() {
  Serial.begin(9600);

  Serial.print("x");
  Serial.print("\t");
  Serial.print("y");
  Serial.print("\t");
  Serial.print("z");
  Serial.print("\t");
  Serial.print("theta");
  Serial.print("\t");
  Serial.print("phi");
  Serial.print("\t");
  Serial.print("psi");
  Serial.println();

  delay(50);
  
  manipulator.forwardKinematics(jointAngles);
  
  Serial.print("Free memory: ");
  Serial.println(freeMemory());
  Serial.println();
}

void loop() {

  manipulator.getPose(pose);

  
  // Angles are printed in degrees.
  // The function calls below shows how easy it is to get the results from the inverse kinematics solution.
  Serial.print(pose[0]);
  Serial.print("\t");
  Serial.print(pose[1]);
  Serial.print("\t");
  Serial.print(pose[2]);
  Serial.print("\t");
  Serial.print(pose[3]);
  Serial.print("\t");
  Serial.print(pose[4]);
  Serial.print("\t");
  Serial.print(pose[5]);
  Serial.println(); Serial.println();



    manipulator.inverseKinematics(desiredPose1[0], 
                                  desiredPose1[1], 
                                  desiredPose1[2], 
                                  desiredPose1[3], 
                                  desiredPose1[4], 
                                  desiredPose1[5], 0, 0.9);


    manipulator.getPose(pose);
  
    // Angles are printed in degrees.
    // The function calls below shows how easy it is to get the results from the inverse kinematics solution.
    Serial.print(pose[0]);
    Serial.print(" --> ");
    Serial.print(desiredPose1[0]);
    Serial.print("\t");
    Serial.print(pose[1]);
    Serial.print(" --> ");
    Serial.print(desiredPose1[1]);
    Serial.print("\t");
    Serial.print(pose[2]);
    Serial.print(" --> ");
    Serial.print(desiredPose1[2]);
    Serial.print("\t");
    Serial.print(pose[3]);
    Serial.print(" --> ");
    Serial.print(desiredPose1[3]);
    Serial.print("\t");
    Serial.print(pose[4]);
    Serial.print(" --> ");
    Serial.print(desiredPose1[4]);
    Serial.print("\t");
    Serial.print(pose[5]);
    Serial.print(" --> ");
    Serial.print(desiredPose1[5]);
    Serial.println();
  
    manipulator.getJointAngles(jointAngles);
    Serial.print(jointAngles[0]);
    Serial.print("\t");
    Serial.print(jointAngles[1]);
    Serial.print("\t");
    Serial.print(jointAngles[2]);
    Serial.print("\t");
    Serial.print(jointAngles[3]);
    Serial.println();



  Serial.print("Free memory: ");
  Serial.println(freeMemory());
  Serial.println();

  Serial.println("============================================");


  while(1);
  
}


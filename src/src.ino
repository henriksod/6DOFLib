
#include "SixDOF.h"
#include "MatrixMath.h"
#include "MemoryFree.h"

double thetas[] = {0, PI/2, -PI/2, 0, PI, 0};
double ds[] = {0.25, 0, 0, 0.2, 0, 0.2};
double alphas[] = {PI/2, 0, -PI/2, PI/2, PI/2, 0};
double as[] = {0, 0.25, 0, 0, 0, 0};

double jointAngles[] = {0, PI/3, PI/2, 0, -PI/4, 0};
double jointAngles2[] = {0, PI/4, 0, 0, PI/3, 0};

double pose[] = {0, 0, 0, 0, 0, 0};

double desiredPose1[] = {-0.4501, 0.15, 0.5777, 2.2923, 1.5234, -0.6599};
double desiredPose2[] = {-0.4485, 0.1225, 0.1147, 2.3934, 2.0215, 0.5019};
double desiredPose3[] = {0, -0.45, 0.05, -1.5708, 1.5708, -2.3562};

SixDOF manipulator(thetas, ds, alphas, as, 6);

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

  do {

    manipulator.inverseKinematics(desiredPose1[0], 
                                  desiredPose1[1], 
                                  desiredPose1[2], 
                                  desiredPose1[3], 
                                  desiredPose1[4], 
                                  desiredPose1[5], 0, 5);

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
    Serial.print("\t");
    Serial.print(jointAngles[4]);
    Serial.print("\t");
    Serial.print(jointAngles[5]);
    Serial.println();

  } while (manipulator.getIKStatus() == SixDOF::FAILED);

  Serial.print("Free memory: ");
  Serial.println(freeMemory());
  Serial.println();

  Serial.println("============================================");
  
  delay(100);
  
  do {

    manipulator.inverseKinematics(desiredPose2[0], 
                                  desiredPose2[1], 
                                  desiredPose2[2], 
                                  desiredPose2[3], 
                                  desiredPose2[4], 
                                  desiredPose2[5], 0.05, 1.8);


    manipulator.getPose(pose);
  
    // Angles are printed in degrees.
    // The function calls below shows how easy it is to get the results from the inverse kinematics solution.
    Serial.print(pose[0]);
    Serial.print(" --> ");
    Serial.print(desiredPose2[0]);
    Serial.print("\t");
    Serial.print(pose[1]);
    Serial.print(" --> ");
    Serial.print(desiredPose2[1]);
    Serial.print("\t");
    Serial.print(pose[2]);
    Serial.print(" --> ");
    Serial.print(desiredPose2[2]);
    Serial.print("\t");
    Serial.print(pose[3]);
    Serial.print(" --> ");
    Serial.print(desiredPose2[3]);
    Serial.print("\t");
    Serial.print(pose[4]);
    Serial.print(" --> ");
    Serial.print(desiredPose2[4]);
    Serial.print("\t");
    Serial.print(pose[5]);
    Serial.print(" --> ");
    Serial.print(desiredPose2[5]);
    Serial.println();
  
    manipulator.getJointAngles(jointAngles);
    Serial.print(jointAngles[0]);
    Serial.print("\t");
    Serial.print(jointAngles[1]);
    Serial.print("\t");
    Serial.print(jointAngles[2]);
    Serial.print("\t");
    Serial.print(jointAngles[3]);
    Serial.print("\t");
    Serial.print(jointAngles[4]);
    Serial.print("\t");
    Serial.print(jointAngles[5]);
    Serial.println();

  } while (manipulator.getIKStatus() == SixDOF::FAILED);

  Serial.print("Free memory: ");
  Serial.println(freeMemory());
  Serial.println();

  Serial.println("============================================");

  delay(100);
  
  do {

    manipulator.inverseKinematics(desiredPose3[0], 
                                  desiredPose3[1], 
                                  desiredPose3[2], 
                                  desiredPose3[3], 
                                  desiredPose3[4], 
                                  desiredPose3[5], 0.05, 1.8);


    manipulator.getPose(pose);
  
    // Angles are printed in degrees.
    // The function calls below shows how easy it is to get the results from the inverse kinematics solution.
    Serial.print(pose[0]);
    Serial.print(" --> ");
    Serial.print(desiredPose3[0]);
    Serial.print("\t");
    Serial.print(pose[1]);
    Serial.print(" --> ");
    Serial.print(desiredPose3[1]);
    Serial.print("\t");
    Serial.print(pose[2]);
    Serial.print(" --> ");
    Serial.print(desiredPose3[2]);
    Serial.print("\t");
    Serial.print(pose[3]);
    Serial.print(" --> ");
    Serial.print(desiredPose3[3]);
    Serial.print("\t");
    Serial.print(pose[4]);
    Serial.print(" --> ");
    Serial.print(desiredPose3[4]);
    Serial.print("\t");
    Serial.print(pose[5]);
    Serial.print(" --> ");
    Serial.print(desiredPose3[5]);
    Serial.println();
  
    manipulator.getJointAngles(jointAngles);
    Serial.print(jointAngles[0]);
    Serial.print("\t");
    Serial.print(jointAngles[1]);
    Serial.print("\t");
    Serial.print(jointAngles[2]);
    Serial.print("\t");
    Serial.print(jointAngles[3]);
    Serial.print("\t");
    Serial.print(jointAngles[4]);
    Serial.print("\t");
    Serial.print(jointAngles[5]);
    Serial.println();

  } while (manipulator.getIKStatus() == SixDOF::FAILED);

  Serial.print("Free memory: ");
  Serial.println(freeMemory());
  Serial.println();

  Serial.println("============================================");

  while(1);
  
}


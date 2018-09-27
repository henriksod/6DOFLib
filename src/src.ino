/********************************************************
 * FABRIK2D 3DOF example
 * Creating the FABRIK object and moving the end effector in a circular motion.
 * You can use whichever unit you want for coordinates, lengths and tolerances as long as you are consistent.
 * Default unit is millimeters.
 ********************************************************/

#include "SixDOF.h"
#include "MatrixMath.h"

float thetas[] = {0, PI/2, -PI/2, 0, PI, 0};
float ds[] = {0.25, 0, 0, 0.2, 0, 0.2};
float alphas[] = {PI/2, 0, -PI/2, PI/2, PI/2, 0};
float as[] = {0, 0.25, 0, 0, 0, 0};

float jointAngles[] = {0, PI/4, 0, 0, PI/4, 0};

float pose[] = {0, 0, 0, 0, 0, 0};

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
}

void loop() {

  Serial.println("Compute const");
  SixDOF manipulator(thetas, ds, alphas, as, 6);
   Serial.println("done");

Serial.println("Compute kin");
  manipulator.forwardKinematics(jointAngles);
Serial.println("done");

  Serial.println("Compute invkin");
  int success = manipulator.inverseKinematics(-0.5182, 0, 0.5682, -PI/2, PI, -PI/2);
  Serial.print(success);
  Serial.println(" done");

  Serial.println();
  
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

  while(1);
}


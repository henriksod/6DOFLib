/********************************************************
 * FABRIK2D 3DOF example
 * Creating the FABRIK object and moving the end effector in a circular motion.
 * You can use whichever unit you want for coordinates, lengths and tolerances as long as you are consistent.
 * Default unit is millimeters.
 ********************************************************/

#include "SixDOF.h"

float thetas[] = {0, PI/2, -PI/2, 0, PI, 0};
float ds[] = {0.25, 0, 0, 0.1, 0, 0.1};
float alphas[] = {PI/2, 0, -PI/2, PI/2, PI/2, 0};
float as[] = {0, 0.25, 0, 0, 0, 0};

SixDOF manipulator(thetas, ds, alphas, as);

float jointAngles[] = {0, 0, 0, 0, 0, 0};

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

  manipulator.forwardKinematics(jointAngles);
}

void loop() {

  manipulator.inverseKinematics(0.4, 0, 0.2, 0, -PI/4, 0);

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
  Serial.print("\t");
  Serial.print(pose[6]);
  Serial.println();

  while(1);
  //delay(50);
}

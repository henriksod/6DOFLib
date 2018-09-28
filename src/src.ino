
#include "SixDOF.h"
#include "MatrixMath.h"

float thetas[] = {0, PI/2, -PI/2, 0, PI, 0};
float ds[] = {0.25, 0, 0, 0.2, 0, 0.2};
float alphas[] = {PI/2, 0, -PI/2, PI/2, PI/2, 0};
float as[] = {0, 0.25, 0, 0, 0, 0};

float jointAngles[] = {0.3, PI/3, PI/2, 0.3, -PI/3, 0};

float pose[] = {0, 0, 0, 0, 0, 0};

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
  
  Serial.println("Compute kin");
  manipulator.forwardKinematics(jointAngles);
  Serial.println("done");
}

void loop() {

  int success = manipulator.inverseKinematics(-0.5165, 0, 0.2018, PI/2, PI, PI/2, 0.001);
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

  /*manipulator.getJointAngles(jointAngles);
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
  Serial.println();*/
}


clear all, close all
l1 = 0.25; l2 = 0.25; l3 = 0.2;
L(1)=Revolute('d', l1, 'a', 0, 'alpha', pi/2);
L(2)=Revolute('d', 0, 'a', l2, 'alpha', 0);
L(3)=Revolute('d', 0, 'a', 0, 'alpha', -pi/2);
L(4)=Revolute('d', l3, 'a', 0, 'alpha', 0);
AngleOffset=[0 pi/2 -pi/2 0];

r=SerialLink(L,'name','6DOF Manipulator Arm','offset',AngleOffset);


%angs = [pi/2 pi/2 pi/2 pi];
angs = [pi/2 pi/4 pi/4 pi];

figure(1);
r.plot(angs);


[T H] = r.fkine(angs);
ax = rotm2axang(T.R)
EUL = tr2eul(T)

T
oc = T.t - l3.*T.R*[0 0 1]'

%%
angs = [pi/2 pi/2 pi/2 pi];
r.plot(angs);
%%
angs = [-1.54	-1.43	-1.73	0.03

];
r.plot(angs);

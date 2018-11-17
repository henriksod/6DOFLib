clear all, close all
l1 = 0.25; l2 = 0.25; l3 = 0.2; l4 = 0.2;
L(1)=Revolute('d', l1, 'a', 0, 'alpha', pi/2);
L(2)=Revolute('d', 0, 'a', l2, 'alpha', 0);
L(3)=Revolute('d', 0, 'a', 0, 'alpha', -pi/2);
L(4)=Revolute('d', l3, 'a', 0, 'alpha', pi/2);
L(5)=Revolute('d', 0, 'a', 0, 'alpha', pi/2);
L(6)=Revolute('d', l4, 'a', 0, 'alpha', 0);
AngleOffset=[0 pi/2 -pi/2 0 pi 0];

r=SerialLink(L,'name','6DOF Manipulator Arm','offset',AngleOffset);

%angs = [0.3 pi/3 pi/2 0.3 -pi/3 0];
%angs = [0 pi/3 pi/2 pi/3 -pi/4 0];
angs = [pi/2 pi/2 pi/2 pi pi/2 pi/4];


%angs = [0 pi/4 0 -pi/3 pi/3 0];

%angs = [0 0 0 0 0.1 0];
%angs = [0 0 0 0 0 pi/2];
figure(1);
r.plot(angs);


[T H] = r.fkine(angs);
ax = rotm2axang(T.R)
EUL = tr2eul(T)

T
oc = T.t - l4.*T.R*[0 0 1]'

%%
angs = [pi/2 pi/2 pi/2 pi pi/2 pi/4];
r.plot(angs);
%%
angs = [1.55	-2.97	-1.98	0.19	0.22	3.74

];
r.plot(angs);

%% DYNMOD Y
% Si vuole giungere all'espressione della matrice regressore.

syms q1 q2 q3 q4 qd1 qd2 qd3 qd4 qdd1 qdd2 qdd3 qdd4
syms ml1 ml2 ml3 m_load Il1 Il2 Il4 Im1 Im2 Im3 Im4 l1 l2 Fm1 Fm2 Fm3 Fm4

a1 = 0.5;
a2 = 0.5;
d0 = 1;
kr1 = 1;
kr2 = 1;
kr3 = 50;
kr4 = 20;

%% calcolo vettori posizione e assi z
p0 = [0 ; 0 ; d0];
p1 = [a1*cos(q1) ; a1*sin(q1) ; d0];
p2 = [a1*cos(q1) + a2*cos(q1+q2) ; a1*sin(q1) + a2*sin(q1+q2) ; d0];
p3 = p2 + [0 ; 0 ; q3]; % C'era un meno
p4 = p3;
pl1 = [l1*cos(q1) ; l1*sin(q1) ; d0];
pl2 = [a1*cos(q1) + l2*cos(q1+q2) ; a1*sin(q1) + l2*sin(q1+q2) ; d0];
pl3 = p3; 
pl4 = p4; 
z0 = [0 0 1]';
z1 = [0 0 1]';
z2 = [0 0 1]';
z3 = [0 0 1]';
z4 = [0 0 1]';

%% calcolo Jacobiani
Jp_l1 = [cross(z0,pl1-p0) zeros(3,3)];
Jp_l2 = [cross(z0,pl2-p0) cross(z1,pl2-p1) zeros(3,2)];
Jp_l3 = [cross(z0,pl3-p0) cross(z1,pl3-p1) z2 zeros(3,1)]; 
Jp_l4 = [cross(z0,pl4-p0) cross(z1,pl4-p1) z2 cross(z3,pl4-p3)]; %m4=mload

Jo_l1 = [z0 zeros(3,3)];
Jo_l2 = [z0 z1 zeros(3,2)];
Jo_l3 = [z0 z1 zeros(3,2)]; %Il3=0
Jo_l4 = [z0 z1 zeros(3,1) z3];


% Le masse dei motori sono trascurati, dunque ci concentriamo sulle parti
% angolari degli jacobiani.

Jo_m1 = Jo_l1;
Jo_m1(:,1) = kr1*z0;
Jo_m2 = Jo_l2;  
Jo_m2(:,2) = kr2*z1;
Jo_m3 = Jo_l3; 
Jo_m3(:,3) = kr3*z2;
Jo_m4 = Jo_l4;  
Jo_m4(:,4) = kr4*z3;


%% calcolo B

B1 = ml1*(Jp_l1'*Jp_l1) + Il1*(Jo_l1'*Jo_l1) + Im1*(Jo_m1'*Jo_m1); 
B2 = ml2*(Jp_l2'*Jp_l2) + Il2*(Jo_l2'*Jo_l2) + Im2*(Jo_m2'*Jo_m2); 
B3 = ml3*(Jp_l3'*Jp_l3) + 0 + Im3*(Jo_m3'*Jo_m3); %Il3=0
B4 = m_load*(Jp_l4'*Jp_l4) + Il4*(Jo_l4'*Jo_l4) + Im4*(Jo_m4'*Jo_m4); 

B = B1 + B2 + B3 + B4;

%% calcolo C

C = [-a1*qd2*sin(q2)*(a2*m_load + a2*ml3 + l2*ml2), -a1*sin(q2)*(qd1 + qd2)*(a2*m_load + a2*ml3 + l2*ml2), 0, 0;
      a1*qd1*sin(q2)*(a2*m_load + a2*ml3 + l2*ml2),  0, 0, 0;
      0, 0, 0, 0;
      0, 0, 0, 0]; 

%% calcolo Fv

Fv = diag([kr1^2*Fm1 kr2^2*Fm2 kr3^2*Fm3 kr4^2*Fm4]);

%% calcolo g

g0 = [0 0 -9.81]';
g1 = - (ml1*g0'*Jp_l1(:,1) + ml2*g0'*Jp_l2(:,1) + ml3*g0'*Jp_l3(:,1) + m_load*g0'*Jp_l4(:,1));
g2 = - (ml1*g0'*Jp_l1(:,2) + ml2*g0'*Jp_l2(:,2) + ml3*g0'*Jp_l3(:,2) + m_load*g0'*Jp_l4(:,2));
g3 = - (ml1*g0'*Jp_l1(:,3) + ml2*g0'*Jp_l2(:,3) + ml3*g0'*Jp_l3(:,3) + m_load*g0'*Jp_l4(:,3));
g4 = - (ml1*g0'*Jp_l1(:,4) + ml2*g0'*Jp_l2(:,4) + ml3*g0'*Jp_l3(:,4) + m_load*g0'*Jp_l4(:,4));
g = [g1 g2 g3 g4]';

%% Calcolo Y

syms qrd1 qrd2 qrd3 qrd4 qrdd1 qrdd2 qrdd3 qrdd4

q=[q1; q2; q3; q4];
qd=[qd1; qd2; qd3; qd4];
qdd=[qdd1; qdd2; qdd3; qdd4];
qrd = [qrd1; qrd2; qrd3; qrd4];
qrdd = [qrdd1; qrdd2; qrdd3; qrdd4];

Tau=B*qrdd+C*qrd+Fv*qrd+g;

Y1_1 = jacobian(Tau, ml1);
Y1_2 = jacobian(Y1_1, l1);
Y1 = jacobian(0.5*Y1_2, l1);
Y2 = jacobian(Tau, ml2);
Y3 = jacobian(Y2, l2);
Y4 = jacobian(0.5*Y3, l2);
Y2 = subs(Y2, l2, 0);
Y3 = subs(Y3, l2, 0);
Y5 = jacobian(Tau, ml3);
Y6 = jacobian(Tau, m_load);
Y7 = jacobian(Tau, Il1);
Y8 = jacobian(Tau, Il2);
Y9 = jacobian(Tau, Il4);
Y10 = jacobian(Tau, Im1);
Y11 = jacobian(Tau, Im2);
Y12 = jacobian(Tau, Im3);
Y13 = jacobian(Tau, Im4);
Y14 = jacobian(Tau, Fm1);
Y15 = jacobian(Tau, Fm2);
Y16 = jacobian(Tau, Fm3);
Y17 = jacobian(Tau, Fm4);
Y = [Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y9 Y10 Y11 Y12 Y13 Y14 Y15 Y16 Y17];

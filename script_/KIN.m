%% Parametri

d0=1;      
a1=0.5;    
a2=0.5;
l1=0.25;
ml1=20;
ml2=20;
ml3=10;
Il1=4;
Il2=4;
Il4=1;
kr1=1;
kr2=1;
kr3=50;
kr4=20;
Im1=0.01;
Im2=0.01;
Im3=0.005;
Im4=0.001;
Fm1=0.00005;
Fm2=0.00005;
Fm3=0.01;
Fm4=0.005;

%% DIRKIN
%% Matrici di trasformazione omogenea
% Per fare in modo di poter successivamente tenere conto della
% tempovarianza di tali matrici le definiremo sulla base di variabili
% simboliche.

syms q1 q2 q3 q4
q=[q1;q2;q3;q4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_B_0=[1 0 0 0 ;
       0 1 0 0 ;
       0 0 1 d0;
       0 0 0 1 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_0_1star=[cos(q1) -sin(q1) 0 0;
           sin(q1)  cos(q1) 0 0;
           0        0       1 0;
           0        0       0 1];
A_1star_1=[1 0 0 a1;
           0 1 0 0 ;
           0 0 1 0 ;
           0 0 0 1 ];
T_0_1=A_0_1star*A_1star_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_1_2star=[cos(q2) -sin(q2) 0 0;
           sin(q2)  cos(q2) 0 0;
           0        0       1 0;
           0        0       0 1];
A_2star_2=[1 0 0 a2;
           0 1 0 0 ;
           0 0 1 0 ;
           0 0 0 1 ];
T_1_2=A_1_2star*A_2star_2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_2_3=[1 0 0 0 ;
       0 1 0 0 ;
       0 0 1 q3;
       0 0 0 1 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_3_4=[cos(q4) -sin(q4) 0 0;
        sin(q4)  cos(q4) 0 0;
        0        0       1 0;
        0        0       0 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_4_E=[1  0  0 0;
       0 -1  0 0;
       0  0 -1 0;
       0  0  0 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_0_4=T_0_1*T_1_2*T_2_3*T_3_4;

T_B_E=T_B_0*T_0_4*T_4_E;

%% INVKYN
%% Calcolo lo Jacobiano unico sfruttando le terne disposte (Vedere testo).

R_0_1=T_0_1(1:3,1:3);
R_1_2=T_1_2(1:3,1:3);
R_2_3=T_2_3(1:3,1:3);
R_3_4=T_3_4(1:3,1:3);

z0=[0;0;1];
z1=R_0_1*z0;
z2=R_0_1*R_1_2*z0;
z3=R_0_1*R_1_2*R_2_3*z0;

p0t=[0;0;0;1];
p1t=T_0_1*p0t;
p2t=T_0_1*T_1_2*p0t;
p3t=T_0_1*T_1_2*T_2_3*p0t;
p4t=T_0_1*T_1_2*T_2_3*T_3_4*p0t;

P0=p0t(1:3);
P1=p1t(1:3);
P2=p2t(1:3);
P3=p3t(1:3);
P4=p4t(1:3);

JP1=cross(z0,P4-P0);
JP2=cross(z1,P4-P1);
JP3=z2;
JP4=cross(z3,P4-P3);

JO1=z0;
JO2=z1;
JO3=[0;0;0];
JO4=z3;

J=[JP1    JP2    JP3    JP4;
   JO1(3) JO2(3) JO3(3) JO4(3)];






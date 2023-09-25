

%% TRAJECTORY PLANNING
% Occorre pianificare la traiettoria dell'EE lungo un percorso definito
% tramite 11 pose dello stesso tale che vi sia almeno un tratto rettilineo,
% uno circolare e 3 punti di via. Per la pianificazione si ricorre a un
% polinomio cubico.

p1=[0.5 0.5 0.75 -pi/2];
p2=[0.5 0 0.5 0];
p3=[0.5 -0.5 0.75 0];
p4=[0 -0.5 0.5 pi/2];
p5=[-0.5 -0.5 0.75 pi/2];
p6=[-0.6 -0.5 0.25 pi];
p7=[-0.6 0 0.25 pi];
p8=[-0.6 0.5 0.25 pi];
p9=[-0.5 0.5 0.75 -pi/2];
p10=[0 0.5 0.5 -pi/2];
p11=p1;
dt=0.4;
ro=0.25;
centro1=[-0.6 -0.25 0.25];
centro2=[-0.6 0.25 0.25];

Poses=[p1; p2; p3; p4; p5; p6; p7; p8; p9; p10; p11];
Ts=0.001;

ti=0;
tf=2;
T=tf*(length(Poses)-1)+5;
T1=0:Ts:T;

CUB_COEFF_POS=[];
CUB_COEFF_PSI=[];
Ti=ti;
Tf=tf;

for i=1:length(Poses)-1
    
    si=0;
    if i==6 || i==7
        sf=0.25*pi;
    else
        sf=norm(Poses(i+1,1:3)-Poses(i,1:3));
    end
    si_dot=0;
    sf_dot=0;

    A=[Ti^3 Ti^2 Ti 1;
        Tf^3 Tf^2 Tf 1;
        3*Ti^2 2*Ti 1 0;
        3*Tf^2 2*Tf 1 0];

    c=[si;sf;si_dot;sf_dot];
    b=A\c;
    CUB_COEFF_POS=[CUB_COEFF_POS; b(1) b(2) b(3) b(4)];
    Ti=Ti+tf;
    Tf=Tf+tf;
end

Ti=ti;
Tf=tf;

for i=1:length(Poses)-1

    si=0;
    sf=norm(Poses(i+1,4)-Poses(i,4));
    si_dot=0;
    sf_dot=0;

    A=[Ti^3 Ti^2 Ti 1;
        Tf^3 Tf^2 Tf 1;
        3*Ti^2 2*Ti 1 0;
        3*Tf^2 2*Tf 1 0];

    c=[si;sf;si_dot;sf_dot];
    b=A\c;
    CUB_COEFF_PSI=[CUB_COEFF_PSI; b(1) b(2) b(3) b(4)];
    Ti=Ti+tf;
    Tf=Tf+tf;
end

S=zeros(10, 2, length(T1));
P_circular=zeros(length(T1),3);
Ti=ti;
Tf=tf;
K=0;

%% Calcolo ascissa curvilinea

for i=1:length(Poses)-1

    a3=CUB_COEFF_POS(i,1);
    a2=CUB_COEFF_POS(i,2);
    a1=CUB_COEFF_POS(i,3);
    a0=CUB_COEFF_POS(i,4);
    a3p=CUB_COEFF_PSI(i,1);
    a2p=CUB_COEFF_PSI(i,2);
    a1p=CUB_COEFF_PSI(i,3);
    a0p=CUB_COEFF_PSI(i,4);
    
    if (i==2) || (i==4) || (i==10)
        K=K+1;
    end
    
    if (i==6)||(i==7)
        for j=1:length(T1)

            if (T1(j)>=Ti-K*dt)&&(T1(j)<=Tf-K*dt)
                S(i,1,j)=a3*(T1(j)+K*dt)^3+a2*(T1(j)+K*dt)^2+a1*(T1(j)+K*dt)+a0;
                S(i,2,j)=a3p*(T1(j)+K*dt)^3+a2p*(T1(j)+K*dt)^2+a1p*(T1(j)+K*dt)+a0p;
            elseif T1(j)>Tf-K*dt
                S(i,1,j)=ro*pi;
                S(i,2,j)=norm(Poses(i+1,4)-Poses(i,4));
            end
        end
    else
        for j=1:length(T1)
            if (T1(j)>=Ti-K*dt)&&(T1(j)<=Tf-K*dt)
                S(i,1,j)=a3*(T1(j)+K*dt)^3+a2*(T1(j)+K*dt)^2+a1*(T1(j)+K*dt)+a0;
                S(i,2,j)=a3p*(T1(j)+K*dt)^3+a2p*(T1(j)+K*dt)^2+a1p*(T1(j)+K*dt)+a0p;
            elseif T1(j)>Tf-K*dt
                S(i,1,j)=norm(Poses(i+1,1:3)-Poses(i,1:3));
                S(i,2,j)=norm(Poses(i+1,4)-Poses(i,4));
            end
        end
    end

    Ti=Ti+tf;
    Tf=Tf+tf;

end

%% Calcolo percorso

 P_e=Poses(1,1:3);
 Psi_e=Poses(1,4);

 for i=1:5
     P_e=P_e+1/norm(Poses(i+1,1:3)-Poses(i,1:3))*S(i,1,:).*(Poses(i+1,1:3)-Poses(i,1:3));
 end

 for i=1:length(T1)
     if T1(i)>10-2*dt && T1(i)<12-2*dt
         P_e(1,:,i) = centro1 + [ro*cos(-S(6,1,i)/ro-pi/2) ro*sin(-S(6,1,i)/ro-pi/2) 0];
     elseif T1(i)>=12-2*dt
         P_e(1,:,i) = Poses(7,1:3);
     end
 end
 for i=1:length(T1)
     if T1(i)>12-2*dt && T1(i)<14-2*dt
         P_e(1,:,i) = centro2 + [ro*cos(S(7,1,i)/ro-pi/2) ro*sin(S(7,1,i)/ro-pi/2) 0];
     elseif T1(i)>=14-2*dt
         P_e(1,:,i) = Poses(8,1:3);
     end
 end

 for i=8:length(Poses)-1
     P_e=P_e+1/norm(Poses(i+1,1:3)-Poses(i,1:3))*S(i,1,:).*(Poses(i+1,1:3)-Poses(i,1:3));
 end

 Psi_e=Psi_e + 1/norm(Poses(2,4)-Poses(1,4))*S(1,2,:).*(Poses(2,4)-Poses(1,4)) + 1/norm(Poses(4,4)-Poses(3,4))*S(3,2,:).*(Poses(4,4)-Poses(3,4))+ 1/norm(Poses(6,4)-Poses(5,4))*S(5,2,:).*(Poses(6,4)-Poses(5,4))+ 1/norm(Poses(9,4)-Poses(8,4))*S(8,2,:).*(Poses(9,4)-Poses(8,4));

 %% Plot

 Pplot=zeros(length(T1),3);
 Psiplot=zeros(length(T1),1);

for i=1:length(T1)
    Pplot(i,1)=P_e(1,1,i);
    Pplot(i,2)=P_e(1,2,i);
    Pplot(i,3)=P_e(1,3,i);
    Psiplot(i,1)=Psi_e(1,1,i);
end

% Plot Percorso geometrico

figure
plot3(Pplot(:,1),Pplot(:,2),Pplot(:,3),'k-');
hold on
plot3(Poses(1,1),Poses(1,2),Poses(1,3),'g.','MarkerSize',20)
for i=2:length(Poses)-1
    plot3(Poses(i,1),Poses(i,2),Poses(i,3),'r.','MarkerSize',20);
end
axis equal
grid on
xlabel('x[m]')
ylabel('y[m]')
zlabel('z[m]')

x_ref=[T1', Pplot, Psiplot];
sim("REF_COMPUTATION.slx");
x_dot_ref=ans.xdot_ref;
x_dotdot_ref=ans.xdotdot_ref;
x_dot=ans.xdot_ref.signals.values;
x_dotdot=ans.xdotdot_ref.signals.values;

% % Plot posizione e psi

% figure
% tlayout = tiledlayout(2, 2);
% nexttile
% plot(T1, x_ref(:,2),'r');
% xlabel('t [s]')
% ylabel('x [m]')
% grid
% nexttile
% plot(T1, x_ref(:, 3),'g');
% xlabel('t [s]')
% ylabel('y [m]')
% grid
% nexttile
% plot(T1, x_ref(:, 4),'b');
% xlabel('t [s]')
% ylabel('z [m]')
% grid
% nexttile
% plot(T1, x_ref(:,5),'m');
% xlabel('t [s]')
% ylabel('psi [rad]')
% grid

% % Plot velocità e psi dot
% 
% figure
% tlayout = tiledlayout(2, 2);
% nexttile
% plot(T1, x_dot(:, 1),'r');
% xlabel('t [s]')
% ylabel('x_d_o_t [m/s]')
% grid
% nexttile
% plot(T1, x_dot(:, 2),'g');
% xlabel('t [s]')
% ylabel('y_d_o_t [m/s]')
% grid
% nexttile
% plot(T1, x_dot(:, 3),'b');
% xlabel('t [s]')
% ylabel('z_d_o_t [m/s]')
% grid
% nexttile
% plot(T1, x_dot(:, 4),'m');
% xlabel('t [s]')
% ylabel('psi_d_o_t [rad/s]')
% grid
% 
% % Plot accelerazione e psi dot dot
% 
% figure
% tlayout = tiledlayout(2, 2);
% nexttile
% plot(T1, x_dotdot(:, 1),'r');
% xlabel('t [s]')
% ylabel('x_d_o_t_d_o_t [m/s]')
% grid
% nexttile
% plot(T1, x_dotdot(:, 2),'g');
% xlabel('t [s]')
% ylabel('y_d_o_t_d_o_t [m/s]')
% grid
% nexttile
% plot(T1, x_dotdot(:, 3),'b');
% xlabel('t [s]')
% ylabel('z_d_o_t_d_o_t [m/s]')
% grid
% nexttile
% plot(T1, x_dotdot(:, 4),'m');
% xlabel('t [s]')
% ylabel('psi_d_o_t_d_o_t [rad/s]')
% grid
% 
% 
% 




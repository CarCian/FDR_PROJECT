
clear
clc

%% CONTROL ALGORITHMS
% Si procede con l'esecuzione degli schemi Simulink realizzanti gli
% algoritmi di controllo robusto, adattativo e ID nello spazio operativo

run("TRAJ.m");

%% Controllo robusto
 
% sim("CTRL_RBST.slx");
% 
% for i=1:length(T1)
% 
%     errxe_rbst(:,i)=ans.ERRXE_RBST.signals.values(:,1,i);
%     xe_rbst(:,i)=ans.XE_RBST.signals.values(:,1,i);
%     errq_rbst(:,i)=ans.ERRQ_RBST.signals.values(i,:);
%     errqd_rbst(:,i)=ans.ERRQD_RBST.signals.values(i,:);
%     u_rbst(:,i)=ans.U_RBST.signals.values(i,:);
% 
% end
% 
% figure
% plot(T1,errxe_rbst(1,:),'r');
% hold on
% plot(T1,errxe_rbst(2,:),'g');
% hold on
% plot(T1,errxe_rbst(3,:),'b');
% hold on
% plot(T1,errxe_rbst(4,:),'m');
% xlabel('t [s]');
% legend('e_x','e_y','e_z','e_p_s_i');
% grid
% 
% figure
% plot3(xe_rbst(1,:),xe_rbst(2,:),xe_rbst(3,:),'r-');
% hold on
% plot3(Poses(1,1),Poses(1,2),Poses(1,3),'g.','MarkerSize',20)
% for i=2:length(Poses)-1
%     plot3(Poses(i,1),Poses(i,2),Poses(i,3),'r.','MarkerSize',20);
% end
% axis equal
% grid on
% xlabel('x[m]')
% ylabel('y[m]')
% zlabel('z[m]')
% 
% figure
% plot(T1,errq_rbst(1,:),'r');
% hold on
% plot(T1,errq_rbst(2,:),'g');
% hold on
% plot(T1,errq_rbst(3,:),'b');
% hold on
% plot(T1,errq_rbst(4,:),'m');
% xlabel('t [s]');
% legend('err q_1','err q_2','err q_3','err q_4');
% grid
% 
% figure
% plot(T1,errqd_rbst(1,:),'r');
% hold on
% plot(T1,errqd_rbst(2,:),'g');
% hold on
% plot(T1,errqd_rbst(3,:),'b');
% hold on
% plot(T1,errqd_rbst(4,:),'m');
% xlabel('t [s]');
% legend('err q_d_o_t _1','err q_d_o_t _2','err q_d_o_t _3','err q_d_o_t _4');
% grid
% 
% figure
% plot(T1,u_rbst(1,:),'r');
% hold on
% plot(T1,u_rbst(2,:),'g');
% hold on
% plot(T1,u_rbst(3,:),'b');
% hold on
% plot(T1,u_rbst(4,:),'m');
% xlabel('t [s]');
% legend('u_1','u_2','u_3','u_4');
% grid

%% Controllo adattativo
% 
% a1 = 0.5;
% a2 = 0.5;
% d0 = 1;
% 
% ml1 = 20;
% ml2 = 20;
% ml3 = 10;
% m_load = 4;
% Il1 = 4;
% Il2 = 4;
% Il4 = 1;
% Im1 = 0.01;
% Im2 = 0.01;
% Im3 = 0.005;
% Im4 = 0.001;
% l1 = 0.25;
% l2 = 0.25;
% kr1 = 1;
% kr2 = 1;
% kr3 = 50;
% kr4 = 20;
% Fm1 =  0.00005;
% Fm2 =  0.00005;
% Fm3 = 0.01;
% Fm4 = 0.005;
% 
% PI0=[ml1*l1^2 ml2 ml2*l2 ml2*l2^2 ml3 0 Il1 Il2 Il4 Im1 Im2 Im3 Im4 Fm1 Fm2 Fm3 Fm4];
% 
% sim("CTRL_ADAPT.slx");
% 
% for i=1:length(T1)
% 
%     errxe_adapt(:,i)=ans.ERRXE_ADAPT.signals.values(:,1,i);
%     xe_adapt(:,i)=ans.XE_ADAPT.signals.values(:,1,i);
%     errq_adapt(:,i)=ans.ERRQ_ADAPT.signals.values(i,:);
%     errqd_adapt(:,i)=ans.ERRQD_ADAPT.signals.values(i,:);
%     u_adapt(:,i)=ans.U_ADAPT.signals.values(i,:);
%     pi_adapt(:,i)=ans.PI_ADAPT.signals.values(:,1,i);
% 
% end
% 
% figure
% plot(T1,errxe_adapt(1,:),'r');
% hold on
% plot(T1,errxe_adapt(2,:),'g');
% hold on
% plot(T1,errxe_adapt(3,:),'b');
% hold on
% plot(T1,errxe_adapt(4,:),'m');
% xlabel('t [s]');
% legend('e_x','e_y','e_z','e_p_s_i');
% grid
% 
% figure
% plot3(xe_adapt(1,:),xe_adapt(2,:),xe_adapt(3,:),'r-');
% hold on
% plot3(Poses(1,1),Poses(1,2),Poses(1,3),'g.','MarkerSize',20)
% for i=2:length(Poses)-1
%     plot3(Poses(i,1),Poses(i,2),Poses(i,3),'r.','MarkerSize',20);
% end
% axis equal
% grid on
% xlabel('x[m]')
% ylabel('y[m]')
% zlabel('z[m]')
% 
% figure
% plot(T1,errq_adapt(1,:),'r');
% hold on
% plot(T1,errq_adapt(2,:),'g');
% hold on
% plot(T1,errq_adapt(3,:),'b');
% hold on
% plot(T1,errq_adapt(4,:),'m');
% xlabel('t [s]');
% legend('err q_1','err q_2','err q_3','err q_4');
% grid
% 
% figure
% plot(T1,errqd_adapt(1,:),'r');
% hold on
% plot(T1,errqd_adapt(2,:),'g');
% hold on
% plot(T1,errqd_adapt(3,:),'b');
% hold on
% plot(T1,errqd_adapt(4,:),'m');
% xlabel('t [s]');
% legend('err q_d_o_t _1','err q_d_o_t _2','err q_d_o_t _3','err q_d_o_t _4');
% grid
% 
% figure
% plot(T1,u_adapt(1,:),'r');
% hold on
% plot(T1,u_adapt(2,:),'g');
% hold on
% plot(T1,u_adapt(3,:),'b');
% hold on
% plot(T1,u_adapt(4,:),'m');
% xlabel('t [s]');
% legend('u_1','u_2','u_3','u_4');
% grid
% 
% figure
% plot(T1,pi_adapt(6,:),'r');
% xlabel('t [s]');
% legend('m_l_o_a_d');
% grid

%% Controllo ID nello spazio operativo

sim('CTRL_IDSO.slx')

for i=1:length(T1)

    errxe_idso(:,i)=ans.ERRXE_IDSO.signals.values(:,1,i);
    errxed_idso(:,i)=ans.ERRXED_IDSO.signals.values(:,1,i);
    xe_idso(:,i)=ans.XE_IDSO.signals.values(:,1,i);
    u_idso(:,i)=ans.U_IDSO.signals.values(i,:);

end

figure
plot(T1,errxe_idso(1,:),'r');
hold on
plot(T1,errxe_idso(2,:),'g');
hold on
plot(T1,errxe_idso(3,:),'b');
hold on
plot(T1,errxe_idso(4,:),'m');
xlabel('t [s]');
legend('e_x','e_y','e_z','e_p_s_i');
grid

figure
plot(T1,errxed_idso(1,:),'r');
hold on
plot(T1,errxed_idso(2,:),'g');
hold on
plot(T1,errxed_idso(3,:),'b');
hold on
plot(T1,errxed_idso(4,:),'m');
xlabel('t [s]');
legend('e_x _d_o_t','e_y _d_o_t','e_z _d_o_t','e_p_s_i _d_o_t');
grid

figure
plot(T1,u_idso(1,:),'r');
hold on
plot(T1,u_idso(2,:),'g');
hold on
plot(T1,u_idso(3,:),'b');
hold on
plot(T1,u_idso(4,:),'m');
xlabel('t [s]');
legend('u_1','u_2','u_3','u_4');
grid

figure
plot3(xe_idso(1,:),xe_idso(2,:),xe_idso(3,:),'r-');
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

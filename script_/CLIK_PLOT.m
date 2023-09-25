%% CLIK ALGORITHMS
% Si procede con l'esecuzione degli schemi Simulink realizzanti i vari
% algoritmi CLIK

run("TRAJ.m");

%% CLIK con Inversa dello Jacobiano

% sim("CLIK_INVJ.slx")
% 
% for i=1:length(T1)
%     err_invj(:,i)=ans.ERR_INVJ.signals.values(:,1,i);
%     q_1(i)=ans.q1_INVJ.signals.values(1,1,i);
%     q_2(i)=ans.q2_INVJ.signals.values(1,1,i);
%     q_3(i)=ans.q3_INVJ.signals.values(1,1,i);
%     q_4(i)=ans.q4_INVJ.signals.values(1,1,i);
% end
% 
% figure
% plot(T1,err_invj(1,:),'r');
% hold on
% plot(T1,err_invj(2,:),'g');
% hold on
% plot(T1,err_invj(3,:),'b');
% hold on
% plot(T1,err_invj(4,:),'m');
% xlabel('t [s]');
% legend('e_x','e_y','e_z','e_p_s_i');
% grid
% 
% figure
% plot(T1, q_1(:) ,'r');
% hold on
% plot(T1, q_2(:) ,'g');
% hold on
% plot(T1, q_3(:) ,'b');
% hold on
% plot(T1, q_4(:) ,'m');
% xlabel('t [s]');
% legend('q_1','q_2','q_3','q_4');
% grid

%% CLIK con Trasposta dello Jacobiano

% sim("CLIK_TRNS.slx")
% 
% for i=1:length(T1)
%     err_trns(:,i)=ans.ERR_TRNS.signals.values(:,1,i);
%     q_1(i)=ans.q1_TRNS.signals.values(1,1,i);
%     q_2(i)=ans.q2_TRNS.signals.values(1,1,i);
%     q_3(i)=ans.q3_TRNS.signals.values(1,1,i);
%     q_4(i)=ans.q4_TRNS.signals.values(1,1,i);
% end
% 
% figure
% plot(T1,err_trns(1,:),'r');
% hold on
% plot(T1,err_trns(2,:),'g');
% hold on
% plot(T1,err_trns(3,:),'b');
% hold on
% plot(T1,err_trns(4,:),'m');
% xlabel('t [s]');
% legend('e_x','e_y','e_z','e_p_s_i');
% grid
% 
% figure
% plot(T1, q_1(:) ,'r');
% hold on
% plot(T1, q_2(:) ,'g');
% hold on
% plot(T1, q_3(:) ,'b');
% hold on
% plot(T1, q_4(:) ,'m');
% xlabel('t [s]');
% legend('q_1','q_2','q_3','q_4');
% grid

%% CLIK con Pseudo-Inversa dello Jacobiano

% sim("CLIK_INVJ.slx")
% 
% W_invj=ans.W_INVJ.signals.values;
% 
% sim("CLIK_PSEUDO.slx")
% 
% for i=1:length(T1)
%     err_pseudo(:,i)=ans.ERR_PSEUDO.signals.values(:,1,i);
%     q_1(i)=ans.q1_PSEUDO.signals.values(1,1,i);
%     q_2(i)=ans.q2_PSEUDO.signals.values(1,1,i);
%     q_3(i)=ans.q3_PSEUDO.signals.values(1,1,i);
%     q_4(i)=ans.q4_PSEUDO.signals.values(1,1,i);
% end
% 
% for i=1:length(T1)
%     xe_pseudo(i,1)=ans.xe_PSEUDO.signals.values(1,1,i);
%     xe_pseudo(i,2)=ans.xe_PSEUDO.signals.values(2,1,i);
%     xe_pseudo(i,3)=ans.xe_PSEUDO.signals.values(3,1,i);
% end
% 
% W_pseudo=ans.W_PSEUDO.signals.values;
% 
% figure
% plot(T1,err_pseudo(2,:),'g');
% hold on
% plot(T1,err_pseudo(3,:),'b');
% hold on
% plot(T1,err_pseudo(4,:),'m');
% xlabel('t [s]');
% legend('e_y','e_z','e_p_s_i');
% grid
% 
% figure
% plot(T1, q_1(:) ,'r');
% hold on
% plot(T1, q_2(:) ,'g');
% hold on
% plot(T1, q_3(:) ,'b');
% hold on
% plot(T1, q_4(:) ,'m');
% xlabel('t [s]');
% legend('q_1','q_2','q_3','q_4');
% grid
% 
% figure
% plot(T1,W_invj(:),'r');
% xlabel('t [s]');
% grid
% 
% figure
% plot(T1,W_pseudo,'bl')
% xlabel('t [s]');
% grid
% 
% figure
% plot3(xe_pseudo(:,1),xe_pseudo(:,2),xe_pseudo(:,3),'k-');
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

%% CLIK del secondo ordine

sim("CLIK_SCND.slx")

for i=1:length(T1)
    err_scnd(:,i)=ans.ERR_SCND.signals.values(:,1,i);
    err_dot_scnd(:,i)=ans.ERR_DOT_SCND.signals.values(:,1,i);
    q_1(i)=ans.q1_SCND.signals.values(1,1,i);
    q_2(i)=ans.q2_SCND.signals.values(1,1,i);
    q_3(i)=ans.q3_SCND.signals.values(1,1,i);
    q_4(i)=ans.q4_SCND.signals.values(1,1,i);
    q_1_dot(i)=ans.q1d_SCND.signals.values(1,1,i);
    q_2_dot(i)=ans.q2d_SCND.signals.values(1,1,i);
    q_3_dot(i)=ans.q3d_SCND.signals.values(1,1,i);
    q_4_dot(i)=ans.q4d_SCND.signals.values(1,1,i);
    q_1_dotdot(i)=ans.q1dd_SCND.signals.values(1,1,i);
    q_2_dotdot(i)=ans.q2dd_SCND.signals.values(1,1,i);
    q_3_dotdot(i)=ans.q3dd_SCND.signals.values(1,1,i);
    q_4_dotdot(i)=ans.q4dd_SCND.signals.values(1,1,i);
end

for i=1:length(T1)
    xe_scnd(i,1)=ans.xe_SCND.signals.values(1,1,i);
    xe_scnd(i,2)=ans.xe_SCND.signals.values(2,1,i);
    xe_scnd(i,3)=ans.xe_SCND.signals.values(3,1,i);
end

figure
plot(T1,err_scnd(1,:),'r');
hold on
plot(T1,err_scnd(2,:),'g');
hold on
plot(T1,err_scnd(3,:),'b');
hold on
plot(T1,err_scnd(4,:),'m');
xlabel('t [s]');
legend('e_x','e_y','e_z','e_p_s_i');
grid

figure
plot(T1,err_dot_scnd(1,:),'r');
hold on
plot(T1,err_dot_scnd(2,:),'g');
hold on
plot(T1,err_dot_scnd(3,:),'b');
hold on
plot(T1,err_dot_scnd(4,:),'m');
xlabel('t [s]');
legend('e_d_o_t _x','e_d_o_t _y','e_d_o_t _z','e_d_o_t _p_s_i');
grid

figure
plot(T1, q_1(:) ,'r');
hold on
plot(T1, q_2(:) ,'g');
hold on
plot(T1, q_3(:) ,'b');
hold on
plot(T1, q_4(:) ,'m');
xlabel('t [s]');
legend('q_1','q_2','q_3','q_4');
grid

figure
plot(T1, q_1_dot(:) ,'r');
hold on
plot(T1, q_2_dot(:) ,'g');
hold on
plot(T1, q_3_dot(:) ,'b');
hold on
plot(T1, q_4_dot(:) ,'m');
xlabel('t [s]');
legend('q_d_o_t _1','q_d_o_t _2','q_d_o_t _3','q_d_o_t _4');
grid

figure
plot(T1, q_1_dotdot(:) ,'r');
hold on
plot(T1, q_2_dotdot(:) ,'g');
hold on
plot(T1, q_3_dotdot(:) ,'b');
hold on
plot(T1, q_4_dotdot(:) ,'m');
xlabel('t [s]');
legend('q_d_o_t_d_o_t _1','q_d_o_t_d_o_t _2','q_d_o_t_d_o_t _3','q_d_o_t_d_o_t _4');
grid

figure
plot3(xe_scnd(:,1),xe_scnd(:,2),xe_scnd(:,3),'bl-');
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
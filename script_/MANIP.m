clc
clear all
close all

run("KIN.m");

%% ANALISI DI MANIPOLABILITA'
%% Nel condurre l'analisi di manipolabilità ci limiteremo a quella della
% struttura portante: questo ci permette di ridurci al caso di un
% manipolatore planare a due bracci. Lavoreremo dunque su uno Jacobiano
% ridotto.

Jr=J(1:2,1:2);

% Configurazioni
p1=[pi/6 -pi/2 0 0];
p2=[pi/2-pi/16 -pi+pi/16 0 0];
p3=[-pi/8 2*pi/8 0 0];
p4=[pi/16 -pi/16 0 0];

P=[p1;p2;p3;p4];

%% Generazione dei grafici

for i=1:size(P,1)

    j=subs(Jr,[q1 q2 q3 q4],P(i,:));
    [V,D]=eig(j*j');
    S=svd(j);

    Cntr=[a1*cos(P(i,1))+a2*cos(P(i,1)+P(i,2)); a1*sin(P(i,1))+a2*sin(P(i,1)+P(i,2))];

    figure
    [xv,yv]=Ellispoid(Cntr,S,V);        % Ellissoide di forza
    plot(xv,yv,'r-');
    hold on
    [xf,yf]=Ellispoid(Cntr,1./S,V);     % Ellissoide di velocità
    plot(xf,yf,'g-');
    hold on
    plot(0, 0, '-o', 'MarkerSize', 10);
    plot([0 a1*cos(P(i,1))],[0 a1*sin(P(i,1))], 'b-');
    plot([a1*cos(P(i,1)) Cntr(1)],[a1*sin(P(i,1)) Cntr(2)], 'b-');
    axis equal
    hold off
    xlabel('x[m]');
    ylabel('y[m]');
    legend('Ellissoide di forza','Ellissoide di velocità','','','');

    w=a1*a2*abs(sin(P(i,2)));
    disp(w)
end



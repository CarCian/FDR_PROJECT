function [x,y] = Ellispoid(C,S,V)

% Plotta un ellissoide sulla base di coordinate del centro, valori
% singolari di J e autovettori di JJ'.

a=S(1);
b=S(2);
t=linspace(0,2*pi);

if (a>=b)
    X=cos(t);
    Y=b/a*sin(t);
elseif (b>a)
    X=0.25*a/b*cos(t);
    Y=0.25*sin(t);
end

ax1=V(:,1);
w=atan2(ax1(2),ax1(1));
x=C(1)+(X*cos(w)-Y*sin(w));
y=C(2)+(X*sin(w)+Y*cos(w));
end
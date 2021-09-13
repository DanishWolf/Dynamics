%% Deflection Equations

P = 9.81*0.5;
l = 0:25:500;
E = 69*10^9;
I = (18.75*3^3)/(12);
Vmax = (5.*l.^3*P)/(24*E*I);
x = 0:.025:.250;
Vmax * 10^3
%v = ((P.*x.^3)/12) - (((P*l^2)/16).*x)
%v * 10^3
%x2 = .275:.025:.500;
%v2 = ((P*l.*x2.^2)/4) - ((P.*x2.^3)/12) - ((3*P*l^2 .*x2)/16) + ((P*l^3)/48)
%v2 * 10^3
%% Deflection Equations

%figure
%plot(x2,v2)

%x_t = 0:25:500;
%v_t = [v v2];

%figure
%plot(x_t,v_t)
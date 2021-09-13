%% Problem Parameters

m = 0.1;                % mass in [kg]
angle = 10;             % initial angle in [deg]
th0 = angle*pi/180;     % initial angled converted to [rad]
thv0 = 0.01;            % initial velocity in [rad/s]
L = 0.5;
Area = 0.7854;
rho = 1.183;
Cdrag = 0.47;
g = 9.81;

%% Analytical Solution

wn = sqrt(g/L);                     % natural frequency [rad/s]
Amp = sqrt(th0^2 +(thv0/wn)^2);     % Amplitude [rad]
phi = atan(thv0/(th0*wn));          % Phase angle [rad]
Period = (2*pi)/wn;

% Parameters

dt = 0.01;          % dt < (1/wn)
t = 10;
Niter = floor(t/dt);

tSim = zeros(1, Niter);
thA = zeros(1, Niter);
dthA = zeros(1, Niter);
ddthA = zeros(1, Niter);

for i = 1:Niter
    tSim(i) = (i-1)*dt;
    tNOW = tSim(i);
    
    thA(i) = Amp*cos(wn*tNOW-phi);
    dthA(i) = -Amp*wn*sin(wn*tNOW-phi);
    ddthA(i) = -Amp*wn^2 *cos(wn*tNOW-phi);
end

plot(tSim, thA);

%% Numerical Solution (NO DRAG)

%Initialise

thN = zeros(1, Niter);        %Angle
thVN = zeros(1, Niter);       %Velocity
thAN = zeros(1, Niter);       %Acceleration

thN(1) = th0-dt*thv0;
thN(2) = th0;

thVN(1) = thv0;
thVN(2) = thv0;

%thNEX = thN(2) + dt*(thv0);

for i = 3:Niter
    
    thNOW = thN(i-1);
    thPRE = thN(i-2);
    
    f = @(thNEX)...
        m*L*(thNEX-2*thNOW+thPRE)/(dt^2)...
        +m*g*sin(thNOW);
    
    thNEX = fzero(f, thNOW);  % solve for 0_N+1
    thN(i) = thNEX;           % Angle (theta)
    thVN(i) = (thNEX - thPRE) ./(2.*dt);  %Velocity
    thAN(i) = (thNEX - 2.*thNOW + thPRE) ./(dt.^2);  %acceleration
    
end

NoDrag_thN = thN;
NoDrag_thVN = thVN;
NoDrag_thAN = thAN;

plot(tSim, NoDrag_thN)

%% Numerical Solution (WITH DRAG)

%Initialise

thN = zeros(1, Niter);        %Angle
thVN = zeros(1, Niter);       %Velocity
thAN = zeros(1, Niter);       %Acceleration

thN(1) = th0-dt*thv0;
thN(2) = th0;

thVN(1) = thv0;
thVN(2) = thv0;

Const = 0.5.*Cdrag.*rho.*Area; %For Drag Force
%thNEX = thN(2) + dt*(thv0);

for i = 3:Niter
    
    thNOW = thN(i-1);
    thPRE = thN(i-2);
    
    f = @(thNEX)...
        m*L*(thNEX-2*thNOW+thPRE)/(dt^2)...
        + Const*(L*((thNEX-thPRE)/(2*dt)))^2*sign(thNEX-thPRE)...
        +m*g*sin(thNOW);
    
    thNEX = fzero(f, thNOW);  % solve for 0_N+1
    thN(i) = thNEX;           % Angle (theta)
    thVN(i) = (thNEX - thPRE) ./(2.*dt);  %Velocity
    thAN(i) = (thNEX - 2.*thNOW + thPRE) ./(dt.^2);  %acceleration
    
end

WithDrag_thN = thN;
WithDrag_thVN = thVN;
WithDrag_thAN = thAN;

plot(tSim, WithDrag_thN)

%% Graph Plots

figure;
plot(tSim, NoDrag_thN); hold on;
plot(tSim, thA);
grid on;

%% Computer Lab Assessment Part 2
% Question 2 Case 1
% Physical Quantities of the System
% The Item
m = 80;              % Total Mass, [kg]
h = 3000;            % Initial Height [m]
v0 = -1;             % Initial Velocity [m/s]
g = -9.81;           % Gravity [m/s^2]

% The Parachute
Amax = 50;           % Maximum Area [m^2]
tact = 5;           % Time Activated [s]
tuf = 60;            % Unfold Time [s]

% Drag
rho0 = 1.183;        % Air Density
Cd = 0.5;            % Drag Coefficient

% Parameters
dt = 0.01;            % dt 
t = 329;
Niter = floor(t/dt);

% Numerical Solution (WITH DRAG)

tN = zeros(1, Niter);
yN = zeros(1, Niter);
dyN = zeros(1, Niter);
ddyN = zeros(1, Niter);

%Initialise
yN(1) = h+dt*v0;
yN(2) = h;

dyN(1) = v0;
dyN(2) = v0;

ddyN(1) = g;
ddyN(2) = g;

tN(1) = dt;
tN(2) = dt + dt;

for i = 3:Niter
    
    tN(i) = i*dt;
    tNow = tN(i);
    
    yNOW = yN(i-1);
    yPRE = yN(i-2);
      
    
    if (tNow <= tact)
        Area = 0;
    elseif (tNow >= tact) && (tNow <= tact + tuf)
        Area = ((tNow - tact)/(tuf))*Amax;
    else
        Area = 50;
    end
    
    f = @(yNEX)...
        m*((yNEX - 2*yNOW + yPRE) /(dt^2))...
        + m*-g...
        + (1/2)*Cd*(rho0*(1 - 10^-4 *yNOW))*Area*((yNEX - yPRE)/(2*dt))^2 * sign((yNEX - yPRE) /(2*dt));
    
    yNEX = fzero(f, yNOW);
    yN(i) = yNEX;           % Height
    dyN(i) = (yNEX - yPRE) /(2*dt);  % Velocity
    yAN(i) = (yNEX - 2*yNOW + yPRE) /(dt^2);  % Acceleration
    
end

plot(tN,dyN)
title('Velocity of the Crate')
xlabel('Time (s)')
ylabel('Velocity (m/s)')


%% Question 2 Case 2

% Physical Quantities of the System
% The Item
m = 100;              % Total Mass, [kg]
h = 3000;            % Initial Height [m]
v0 = -1;             % Initial Velocity [m/s]
g = -9.81;           % Gravity [m/s^2]

% The Parachute
Amax = 50;           % Maximum Area [m^2]
tact = 10;           % Time Activated [s]
tuf = 30;            % Unfold Time [s]

% Drag
rho0 = 1.183;        % Air Density
Cd = 0.5;            % Drag Coefficient

% Parameters
dt = 0.01;            % dt 
t = 257;
Niter = floor(t/dt);

% Numerical Solution (WITH DRAG)

tN = zeros(1, Niter);
yN = zeros(1, Niter);
dyN = zeros(1, Niter);
ddyN = zeros(1, Niter);

%Initialise
yN(1) = h+dt*v0;
yN(2) = h;

dyN(1) = v0;
dyN(2) = v0;

ddyN(1) = g;
ddyN(2) = g;

tN(1) = dt;
tN(2) = dt + dt;

for i = 3:Niter
    
    tN(i) = i*dt;
    tNow = tN(i);
    
    yNOW = yN(i-1);
    yPRE = yN(i-2);
      
    
    if (tNow <= 10)
        Area = 0;
    elseif (tNow >= tact) && (tNow <= tact + tuf)
        Area = ((tNow - tact)/(tuf))*Amax;
    else
        Area = 50;
    end
    
    f = @(yNEX)...
        m*((yNEX - 2*yNOW + yPRE) /(dt^2))...
        + m*-g...
        + (1/2)*Cd*(rho0*(1 - 10^-4 *yNOW))*Area*((yNEX - yPRE)/(2*dt))^2 * sign((yNEX - yPRE) /(2*dt));
    
    yNEX = fzero(f, yNOW);
    yN(i) = yNEX;           % Height
    dyN(i) = (yNEX - yPRE) /(2*dt);  % Velocity
    yAN(i) = (yNEX - 2*yNOW + yPRE) /(dt^2);  % Acceleration
    
end

plot(tN,dyN)
title('Velocity of the Crate')
xlabel('Time (s)')
ylabel('Velocity (m/s)')



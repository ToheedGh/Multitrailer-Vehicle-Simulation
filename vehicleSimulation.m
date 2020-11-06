% A simulation example of the generated vehicle equations of motion using
% the residual function 'equationsOfMotion.m'.

% Reference publication:
% T. Ghandriz, B. Jacobson, P. Nilsson, L. Laine, and N. Fröjd, “Computationally efficient
% nonlinear one- and two-track models for multitrailer road vehicles,”
% IEEE Access, vol. -, pp. –, 2020.

% Written by: Toheed Ghandriz, April 2020

clear all
close all
clc

parameters; % reading parameters

Animate = 'Yes'; 

% solver settings: 
y0 = [0; 0; 0;0;0;0;param.vx0;0;0;0;0;0];  % initial values of the states
yp0 = [0; 0; 0;0;0;0;0;0;0;0;0;0];         % initial values of the state derivatives
fixedInit = [0; 0; 0;0;0;0;1;0;0;0;0;0].'; % fixed initial values of the states
options = odeset('RelTol',1e-6,'AbsTol',zeros(1,length(y0))+1e-4);
[y0,yp0] = decic(@(t,y,yp) equationsOfMotion(t,y,yp,param),0,y0,fixedInit,yp0,[],options); % fininding consistent initial conditions
tic
[t,y] = ode15i(@(t,y,yp) equationsOfMotion(t,y,yp,param),param.tspan,y0,yp0,options);      % integrating the implicit ODE
toc

% states:
vx = y(:,7);
vy = y(:,8);
phi1=y(:,3);
phi1dot=y(:,9);
thet1=y(:,4);
thet1dot = y(:,10);
thet2=y(:,5);
thet2dot = y(:,11);
thet3=y(:,6);
thet3dot = y(:,12);

X0 = 0; % initial global X-position
Y0 = 0; % initial global Y-position
% transforming from vehicle local to global coordinates
Vx = vx.*cos(y(:,3))-vy.*sin(y(:,3));
Vy = vx.*sin(y(:,3))+vy.*cos(y(:,3));
% integrating the global speeds to find the vehicle trajectory:
X = [X0;cumsum((Vx(1:end-1)+Vx(2:end))/2.*diff(t))+X0];
Y = [Y0;cumsum((Vy(1:end-1)+Vy(2:end))/2.*diff(t))+Y0];


YRu1 = phi1dot;                            % yaw rate unit 1
YRu2 = phi1dot-thet1dot;                   % yaw rate unit 2
YRu4 = phi1dot-thet1dot-thet2dot-thet3dot; % yaw rate unit 4


%% Plotting the trajectory
fig =figure(1);
plot(X,Y)
xlabel('X (m)')
ylabel('Y (m)')
title('Vehicle trajectory')
grid on
hold on

%% Plotting speed and road wheel steering angle
fig =figure(2);
yyaxis right
param.del1_1(1) = 0;
plot(t,param.del1Ref(1:length(t)),'-.')
xlabel('Time (s)')
ylabel('Road wheel steering angle, \delta_1_1 (rad)')
grid on
set(gca,'YColor',[0 0 0]);
xlim([0 param.tEnd])

yyaxis left
plot(t,param.vxRef(1:length(t))*3.6,'-')
hold on
plot(t,vx*3.6,'k--')
xlabel('Time (s)');
ylabel('Longitudinal velocity (km/h)');
legend('Velocity ref','Velocity followed','Steering angle input','location','southeast')
ax = gca;
ax.FontSize = 9;
set(gca,'YColor',[0 0 0]);
xlim([0 param.tEnd])
ylim([0 max(param.vxRef)*3.6*1.3])

fig.Position = [488.0000  300  560.0000  330];


%% Plotting yaw rates
fig = figure(3);
a=subplot('Position',[.1 .70 .85 .25]);
param.testTyreRelaxation = 0;
plot(t,YRu1,'LineWidth',1)
title('6-axle 1-track, nonlinear vehicle')
ax = gca;
ax.FontSize = 11;
ylabel('YR 1st unit (rad/s)')
grid on
xlim([0 param.tEnd])
hold on

b = subplot('Position',[.1 .40 .85 .25]);
plot(t,YRu2, 'LineWidth',1)
ylabel('YR 2nd unit (rad/s)')
xlim([0 param.tEnd])
ax = gca;
ax.FontSize = 11;
grid on
hold on

c = subplot('Position',[.1 .10 .85 .25]);
plot(t,YRu4, 'LineWidth',1)
xlabel('Time (s)')
ylabel('YR 4th unit (rad/s)')
xlim([0 param.tEnd])
grid on
ax = gca;
ax.FontSize = 11;
hold off
fig.Position = 1.0e+02 *[3.802   2.25   6.678   5.37];

%% Animation 
if strcmp(Animate, 'Yes')
    animateSimulation;
end
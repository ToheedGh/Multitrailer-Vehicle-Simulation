% Definition of the parameters (and the default inputs) of the vehicle equation of motion
% These parameters are valid for a single-track 4-unit 6-axle vehicle.
% Reference publication:
% T. Ghandriz, B. Jacobson, P. Nilsson, L. Laine, and N. Fröjd, “Computationally efficient
% nonlinear one- and two-track models for multitrailer road vehicles,”
% IEEE Access, vol. -, pp. –, 2020.

% Written by: Toheed Ghandriz, April 2020


nUnits = 4;           % number of vehicle units
nAxles = 2;           % number of axles of the first vehicle unit 

%% Default inputs
param.Fxw = zeros(nUnits,nAxles);   % (N) axle forces 
param.del = zeros(nUnits,nAxles);   % (rad) road-wheels steering angles 

%% Vehicle parameters

param.m = [10250; 22000; 3700; 24300];                       % (kg) vehicle units masses 
param.J = 1.0e+05 *[0.467116; 1.629069; 0.070802; 3.089547]; % (kg m^2) vehicle units moments of inertia 
param.Fz0A= 1.0e+05 *[ 0.68867    0.95829
                       1.55474         0
                       0.53410    0.53410
                       1.64063         0];    % (N) axle static vertical load on a flat ground

param.xC = [ 0        -2.2008
             5.1313   -5.3187
             3.8999   -0.4001
             5.4222         0];               % (m) vehicle unit coupling local x position w.r.t. unit's COG (first column corresponds to the position of the front coupling, and the second column gives the position of the rear coupling).

param.mu = ones(nUnits,nAxles);               % road friction at each units axle (could potentially be a function of the distance traveled.)
param.xA = [ 1.4992   -2.5858
            -2.5687         0
             0.5499   -1.4501
            -2.4778         0];               % (m) local x position of the axle w.r.t. unit's COG 

param.nW = [2     8
            6     0
            2     2
            6     0];                             % Number of wheels per axle

param.Fz0wNom = 25000 + zeros(nUnits,nAxles);     % (N) wheel nominal vertical load (could be different for different axles.)

param.e = 1;                                      % friction ellipse scale factor (if e is 1 then the ellipse becomes a circle)
param.fr = 0.008;                                 % rolling resistance coefficient

% parameters needed for lateral load transfer:
param.h = [0.63; 1.82; 0; 1.36];                  % (m) units COG height
param.tW = [2.09   2.09
            2.05     0
            2.05   2.05
            2.10     0];                          % (m) axles track width 
param.hRCA = [0.826   1.660
              1.040      0
              1.040  1.040
              1.040      0];                      % (m) axles roll center height 

param.crollA = 1.0e+06*[ 0.468136   0.939976
                         4.398107          0
                         1.468055  1.468055
                         4.409051          0 ];   % (N m/rad) axle rolling stiffness 

% Nonlinear tyre model tuning parameters, different for steerable wheels
% and unsteerable (normal) wheels (tuned for the case of no lateral load transfer):
param.uyg_steerWh = -0.168122; 
param.CCy0_steerWh = 5.33168;
param.uyg_normalWh = -0.10;
param.CCy0_normalWh = 12.38360;

% Linear tyre model parameters:
param.C  = 1.0e+06 * [ 0.167665   0.108188
                       0.474479          0
                       0.494192       1.50
                       0.142017          0];   % (N) axles cornering stiffness 

%% Environment parameters
param.g = 9.810;                  % (m/s^2) Gravitational acceleration
param.AfcdRhoa = 9.9840;          % (N/(m/s)^2) air resistance coefficients
param.pitchAng = zeros(nUnits,1); % (rad) constant road pitch angles at cog of each unit (positive downhill)(could potentially be a function of the distance traveled.)
param.bankAng = zeros(nUnits,1);  % (rad) constant Road banking angles around local x-axis, positive downhill at the left side (could potentially be a function of the distance traveled.)

%% Maneuver parameters

param.sampleRate = 100;                         % (Hz) sample freauency 
param.tEnd = 60;                                % (s) simulation end time 
nSamp = param.tEnd*param.sampleRate+1;          % total number of samples
param.tspan = 0:1/param.sampleRate:param.tEnd;  % (s) time vector
param.vxRef = 30/3.6 * ones(nSamp,1);           % (m/s) desired reference longitudinal speed
t = param.tspan;
steerFreq = 1/20;                               % (Hz) sin steering frequency 
steerAmp = 0.2;                                % (Hz) sin steering frequency
param.del1Ref = steerAmp*sin(2*pi*t*steerFreq)'; % (rad) steering angle reference 
% param.del1Ref = steerAmp*ones(size(param.del1Ref))'; % (rad) steering angle reference

param.vx0 = param.vxRef(1);                     % (m/s) initial longitudinal speed (the initial value of all the other states are zero.)


%% Animation parameters

param.len = [5.3,12,4,13.6];      % Units lengths
param.wt = [2.1,2.1,2.1,2.1];     % Units track width
param.cog2R = [4,6.17,0.95,6.78]; % COG distance from rear edge


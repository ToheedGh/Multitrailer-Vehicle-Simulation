% This code generates all the differential equations describing the dynamic motion of a
% multitrailer vehicle. 
% This code was generated using MATLAB R2019b and does not work in earlier
% versions.

% Reference publication:
% T. Ghandriz, B. Jacobson, P. Nilsson, L. Laine, and N. Fröjd, “Computationally efficient
% nonlinear one- and two-track models for multitrailer road vehicles,”
% IEEE Access, vol. -, pp. –, 2020.

% Written by: Toheed Ghandriz, April 2020

close all
clear all
clc

 
%% Options

singleTrack = 1;                           % 0 or false: tow-track, 1 or true: single-track ('singleTrack = 1': default vehicle)
linearTyre = 1;                            % 0 or false: nonlinear tyre, 1 or true: linear tyre
includeLateralLoadTransfer = 0; 
extractForceElementsFromEqs = 0;           % set as true for a simulation purpose in order to reduce the size of the equations (not necessarily faster simulations). Set as false with the purpose of Jacobian evaluation and nonlinear model predictive control
includeCombinedSlip = 1;                      
includeRoadGradeForce = 1;
includeAirResistance = 1;
includeRollingResistance = 0;

% % default vehicle, a four-unit vehicle with six axles (the first and forth axles are steerable, and single-track)
ua = [1,1;1,0;1,1;1,0];               % (Units and Axles binary matrix) each row has to have the same number of elements, write 0 if an axle does not exist. First unit should have at least 2 axels
sa = [1,0;0,0;1,0;0,0];               % steerable axles

% % Some examples of vehicles: (For the below examples to work, files 'parameters.m' and 'vehicleSimulation.m' need to be updated.)

% % A single unit vehicle with two axles (the first axle is steerable)
% ua = [1,1];                         % (Units and Axles binary matrix) each row has to have the same number of elements, write 0 if an axle does not exist. First unit should have at least 2 axels
% sa = [1,0];                         % steerable axles

% % A two-unit vehicle with three axles (the first axle is steerable)
% ua = [1,1;1,0];                     % (Units and Axles binary matrix) each row has to have the same number of elements, write 0 if an axle does not exist. First unit should have at least 2 axels
% sa = [1,0;0,0];                     % steerable axles

% % A four-unit vehicle with 11 axles (the first and seventh axles are steerable)
% ua  = [1,1,1;1,1,1;1,1,0;1,1,1];    % (Units and Axles binary matrix) each row has to have the same number of elements, write 0 if an axle does not exist. First unit should have at least 2 axels
% sa = [1,0,0;0,0,0;1,0,0;0,0,0];     % steerable axles 


%% Defining symbolic variables and parameters
disp('Generating equations ... ')
disp('It takes some time (less than a minute) ...')

nUnits = size(ua,1);
nUnitAxles = size(ua,2);
syms e ...                              % friction circle scale factor (if e is not 1 then the circle changes to an ellipse)
uyg_steerWh CCy0_steerWh uyg_normalWh...
CCy0_normalWh g AfcdRhoa fr

                                               

% The parameters below can be time dependent or not. Their time derivative
% do not show up in equatoins, so they can have different values in each time step. 

xA = sym('xA',[nUnits nUnitAxles]);       % local x position of the axle w.r.t. unit COG
tW = sym('tW',[nUnits nUnitAxles]);       % track Width
xC = sym('xC',[nUnits 2]);                % unit coupling local x position w.r.t. unit COG
m = sym('m',[nUnits 1]);                  % mass of units   
J = sym('J',[nUnits 1]);                  % inertia of units 
h = sym('h',[nUnits 1]);                  % unit COG height
nW = sym('nW',[nUnits nUnitAxles]);       % number of wheels per axle
hRCA= sym('hRCA',[nUnits nUnitAxles]);    % roll center height for all axles
crollA = sym('crollA',[nUnits nUnitAxles]); % Axle rolling stiffness 
Fz0A = sym('Fz0A',[nUnits nUnitAxles]);   % Axle static vertical load
Fz0wNom = sym('Fz0wNom',[nUnits nUnitAxles]);   % wheel nominal vertical load
pitchAng = sym('pitchAng',[nUnits 1]);    % road pitch angle at cog of each unit (positive downhill)
bankAng = sym('bankAng',[nUnits 1]);      % road banking angle around local x-axis, positive downhill at the left (driver) side

FzlLt = zeros(size(ua)); % lateral load transfer of the left-side axle vertical load (zeros if lateral load transfer is not applied)
FzrLt = zeros(size(ua)); % lateral load transfer of the right-side axle vertical load (zeros if lateral load transfer is not applied)
FzlLt = num2cell(FzlLt);
FzrLt = num2cell(FzrLt);

% states
syms X(t) Y(t) vx(t) vy(t) derVx_t(t) ...
    derVy_t(t) phi1(t) derPhi_t1(t) ...
    derderPhi_t1(t)                       % states of the first unit, (X,Y): global position
syms thet(t) [1 nUnits-1]                 % articulation angles
thet=thet(t);
syms derThet_t(t)  [1 nUnits-1]           % time derivative of articulation angles
derThet_t = derThet_t(t);
syms derderThet_t(t)  [1 nUnits-1]        % time derivative of articulation angles rates
derderThet_t = derderThet_t(t);

if singleTrack
    syms del(t) [nUnits nUnitAxles]       % defining steering angles one for each wheel
    syms C [nUnits nUnitAxles]            % A single wheel cornering stiffness
    syms Fxw(t) [nUnits nUnitAxles]       % side or axle longitudinal force in wheel frame
    syms mu [nUnits nUnitAxles]           % road friction coefficient (can be different for different axles and sides depending on road position)
else
    syms del(t) [nUnits 2*nUnitAxles] 
    syms C [nUnits 2*nUnitAxles]
    syms Fxw(t) [nUnits 2*nUnitAxles]     % side or axle longitudinal force in wheel frame
    syms mu [nUnits 2*nUnitAxles]
end
del = del(t);
% Input forces
Fxw = Fxw(t);

%% Articulation angles and their derivatives
phi{1} = phi1;
derPhi_t{1} = derPhi_t1;
if nUnits>1 
    for i=2:nUnits                     % loop over units
        phi{i} = phi{i-1} - thet(i-1); % According to sign convention defined in ISO-8855
        derPhi_t{i} = derPhi_t{i-1} - derThet_t(i-1);
    end
end


%% Global COG velocities:

derX_t = cell(nUnits,1); 
derY_t = cell(nUnits,1); 
derX_t{1} = cos(phi{1})*vx - sin(phi{1})*vy; % converting local velocities of the first unit to global velocities
derY_t{1} = sin(phi{1})*vx + cos(phi{1})*vy;

if nUnits>1                                  % Global velocities of COG of other units (could be also found by taking the time derivative of COG position)
    for i=2:nUnits 
        derX_t{i} = derX_t{i-1} - xC(i-1,2)*sin(phi{i-1})*derPhi_t{i-1} + xC(i,1)*sin(phi{i})*derPhi_t{i};
        derY_t{i} = derY_t{i-1} + xC(i-1,2)*cos(phi{i-1})*derPhi_t{i-1} - xC(i,1)*cos(phi{i})*derPhi_t{i};
    end
end
                          


%% Kinetic energy
T = 0;
for i=1:nUnits  
    T = T + 0.5*(J(i)*derPhi_t{i}^2 + m(i)*(derX_t{i}^2+derY_t{i}^2));
end

% We know that:
% dT/dX = 0
% dT/dY = 0

A = functionalDerivative(T,vx);
B = functionalDerivative(T,vy);

derA_t = diff(A,t);
derB_t = diff(B,t);

derA_t = subs(derA_t,diff(vx(t),t),derVx_t(t));
derA_t = subs(derA_t,diff(vy(t),t),derVy_t(t));
derA_t = subs(derA_t,diff(phi1(t),t),derPhi_t1(t));
derA_t = subs(derA_t,diff(derPhi_t1(t),t),derderPhi_t1(t));

derB_t = subs(derB_t,diff(vx(t),t),derVx_t(t));
derB_t = subs(derB_t,diff(vy(t),t),derVy_t(t));
derB_t = subs(derB_t,diff(phi1(t),t),derPhi_t1(t));
derB_t = subs(derB_t,diff(derPhi_t1(t),t),derderPhi_t1(t));

if nUnits>1
    for i=1:nUnits-1 
        derA_t = subs(derA_t,diff(thet(i),t),derThet_t(i));
        derA_t = subs(derA_t,diff(derThet_t(i),t),derderThet_t(i));
        derB_t = subs(derB_t,diff(thet(i),t),derThet_t(i));
        derB_t = subs(derB_t,diff(derThet_t(i),t),derderThet_t(i));
    end
end

A = simplify(A);
B = simplify(B);
derA_t = simplify(derA_t);
derB_t = simplify(derB_t);

Eq_A=['A = ',char(A),';']; Eq_A = erase(Eq_A,"(t)"); % converting to string and removing '(t)'
Eq_B=['B = ',char(B),';']; Eq_B = erase(Eq_B,"(t)");
Eq_derA_t=['derA_t = ',char(derA_t),';']; Eq_derA_t = erase(Eq_derA_t,"(t)");
Eq_derB_t=['derB_t = ',char(derB_t),';']; Eq_derB_t = erase(Eq_derB_t,"(t)");

%% Global positions of forces: 

pair = [X + xA(1,1)*cos(phi1); Y + xA(1,1)*sin(phi1)];

% Position of COGs:
pcog = cell(nUnits,1);
pcog{1}=[X;Y];
if nUnits>1
    for i = 2:nUnits
       pcog{i} =  pcog{i-1}+ [xC(i-1,2)*cos(phi{i-1})- xC(i,1)*cos(phi{i});...
                              xC(i-1,2)*sin(phi{i-1})- xC(i,1)*sin(phi{i})];
    end
end  


if singleTrack
    pW = cell(nUnits,nUnitAxles); % wheel positions
    for i = 1:nUnits
        for j=1:nUnitAxles 
            if ua(i,j)>0          % if axle exists, otherwise its position is empty
                pW{i,j} =  pcog{i}+ [xA(i,j)*cos(phi{i});...
                                     xA(i,j)*sin(phi{i})];
            end
        end
    end
else
    pW = cell(nUnits,nUnitAxles*2); % wheel positions (the first half of columns contains the position of the right wheels and the second half contains the position of the left wheels.)
    for i = 1:nUnits
        for j=1:nUnitAxles
            if ua(i,j)>0 % if axle exists
                pW{i,j} =             pcog{i}+ [xA(i,j)*cos(phi{i}) + tW(i,j)/2*sin(phi{i});...  % right wheels
                                                xA(i,j)*sin(phi{i}) - tW(i,j)/2*cos(phi{i})];
                pW{i,j+nUnitAxles} =  pcog{i}+ [xA(i,j)*cos(phi{i}) - tW(i,j)/2*sin(phi{i});...  % left wheels
                                                xA(i,j)*sin(phi{i}) + tW(i,j)/2*cos(phi{i})];
            end
        end
    end
end
    

%% Global velocities at the position of COGs and wheels:

if singleTrack
    derPW_t = cell(nUnits,nUnitAxles); % wheel positions
    for i = 1:nUnits
        for j=1:nUnitAxles
            if ua(i,j)>0 % if axle exists
                derPW_t{i,j} = diff(pW{i,j},t);
                derPW_t{i,j} = subs(derPW_t{i,j},diff(X(t),t),derX_t{1});
                derPW_t{i,j} = subs(derPW_t{i,j},diff(Y(t),t),derY_t{1});
                derPW_t{i,j} = subs(derPW_t{i,j},diff(phi1(t),t),derPhi_t1(t));
                if nUnits>1
                    for k = 1:nUnits-1
                        derPW_t{i,j} = subs(derPW_t{i,j},diff(thet(k),t),derThet_t(k));
                    end
                end
            end
        end
    end
else % if two-track
    derPW_t = cell(nUnits,2*nUnitAxles); % wheel positions
    for i = 1:nUnits
        for j=1:nUnitAxles
            if ua(i,j)>0 % if axle exists
                derPW_t{i,j} = diff(pW{i,j},t);
                derPW_t{i,j} = subs(derPW_t{i,j},diff(X(t),t),derX_t{1});
                derPW_t{i,j} = subs(derPW_t{i,j},diff(Y(t),t),derY_t{1});
                derPW_t{i,j} = subs(derPW_t{i,j},diff(phi1(t),t),derPhi_t1(t));
                if nUnits>1
                    for k = 1:nUnits-1
                        derPW_t{i,j} = subs(derPW_t{i,j},diff(thet(k),t),derThet_t(k));
                    end
                end
                
                derPW_t{i,j+nUnitAxles} = diff(pW{i,j+nUnitAxles},t);
                derPW_t{i,j+nUnitAxles} = subs(derPW_t{i,j+nUnitAxles},diff(X(t),t),derX_t{1});
                derPW_t{i,j+nUnitAxles} = subs(derPW_t{i,j+nUnitAxles},diff(Y(t),t),derY_t{1});
                derPW_t{i,j+nUnitAxles} = subs(derPW_t{i,j+nUnitAxles},diff(phi1(t),t),derPhi_t1(t));
                if nUnits>1
                    for k = 1:nUnits-1
                        derPW_t{i,j+nUnitAxles} = subs(derPW_t{i,j+nUnitAxles},diff(thet(k),t),derThet_t(k));
                    end
                end
                
            end
        end
    end
end

          

%% Velocities at wheel positions in wheel frames:

if singleTrack
    derPW_tWCoor = cell(nUnits,nUnitAxles); % wheel positions
    for i = 1:nUnits
        for j=1:nUnitAxles
            if ua(i,j)>0 % if axle exists
                if sa(i,j)>0  % if axle steerable
                    wheelAngle = phi{i} + del(i,j); 
                else
                    wheelAngle = phi{i};
                end
                derPW_tWCoor{i,j} = [cos(wheelAngle)  , sin(wheelAngle);
                                    -sin(wheelAngle) , cos(wheelAngle)] * derPW_t{i,j};
                derPW_tWCoor{i,j} = derPW_tWCoor{i,j}(t);
            end 
        end
    end
else
    derPW_tWCoor = cell(nUnits,2*nUnitAxles); % wheel positions
    for i = 1:nUnits
        for j=1:nUnitAxles
            if ua(i,j)>0      % if axle exists
                if sa(i,j)>0  % if axle steerable
                    wheelAngle = phi{i} + del(i,j); % right wheel angle (each wheel can have its own angle, left or right can have different angles)
                else
                    wheelAngle = phi{i};
                end
                derPW_tWCoor{i,j} = [cos(wheelAngle)  , sin(wheelAngle);
                                    -sin(wheelAngle) , cos(wheelAngle)] * derPW_t{i,j};
                derPW_tWCoor{i,j} = derPW_tWCoor{i,j}(t);
                                
                if sa(i,j)>0  % if axle steerable
                    wheelAngle = phi{i} + del(i,j+nUnitAxles); % each wheel can have its own angle, left or right can have different angles
                else
                    wheelAngle = phi{i};
                end
                derPW_tWCoor{i,j+nUnitAxles} = [cos(wheelAngle)  , sin(wheelAngle);
                                               -sin(wheelAngle) , cos(wheelAngle)] * derPW_t{i,j+nUnitAxles};
                derPW_tWCoor{i,j+nUnitAxles}= derPW_tWCoor{i,j+nUnitAxles}(t);
            end 
        end
    end
end

%% Calculating lateral acc of cog of the units in global coord (needed for lateral load transfer or post processing of the simulation results. They are not printed in the output file.)

derderY_t= cell(nUnits,1);
for i = 1:nUnits
    derderY_t{i} = diff(derY_t{i},t);
    
    derderY_t{i} = subs(derderY_t{i},diff(vx(t),t),derVx_t(t));
    derderY_t{i} = subs(derderY_t{i},diff(vy(t),t),derVy_t(t));
    derderY_t{i} = subs(derderY_t{i},diff(phi1(t),t),derPhi_t1(t));
    derderY_t{i} = subs(derderY_t{i},diff(derPhi_t1(t),t),derderPhi_t1(t));
    
    if nUnits>1
        for k = 1:nUnits-1
            derderY_t{i} = subs(derderY_t{i},diff(thet(k),t),derThet_t(k));
            derderY_t{i} = subs(derderY_t{i},diff(derThet_t(k),t),derderThet_t(k));
        end
    end
end

derderX_t= cell(nUnits,1);
for i = 1:nUnits
    derderX_t{i} = diff(derX_t{i},t);
    
    derderX_t{i} = subs(derderX_t{i},diff(vx(t),t),derVx_t(t));
    derderX_t{i} = subs(derderX_t{i},diff(vy(t),t),derVy_t(t));
    derderX_t{i} = subs(derderX_t{i},diff(phi1(t),t),derPhi_t1(t));
    derderX_t{i} = subs(derderX_t{i},diff(derPhi_t1(t),t),derderPhi_t1(t));
    
    if nUnits>1
        for k = 1:nUnits-1
            derderX_t{i} = subs(derderX_t{i},diff(thet(k),t),derThet_t(k));
            derderX_t{i} = subs(derderX_t{i},diff(derThet_t(k),t),derderThet_t(k));
        end
    end
end

derv_t= cell(nUnits,1);  % acceleration of the unit COG relative to the unit itself
dervy_t= cell(nUnits,1); % unit COG local lateral acceleration
for i=1:nUnits
    derv_t{i} = [cos(phi{i}) , sin(phi{i});
                 -sin(phi{i}) , cos(phi{i})] * [derderX_t{i};derderY_t{i}];
    derv_t{i} = derv_t{i}(t);
    dervy_t{i} = simplify(derv_t{i}(2)); % 'simplify' can take some time to be executed!!!
end


%% Laterl load transfer (LLT), roll-free couplings

if includeLateralLoadTransfer
    
    nAxUnit = sum(ua.');        % number of active axles per unit
    % LLT for the first unit
    lr = -sum(xA(1,2:end))/(nAxUnit(1)-1); % distance to rear axle group of the first unit. (First unit should have at least 2 axels)
    lf = xA(1,1);
    L = lf+lr;
    hrg1 = sum(hRCA(1,2:end))/(nAxUnit(1)-1); % group-axle roll center height
    crg1 = sum(crollA(1,2:end));              % group-axle roll stiffness
    dh = h(1) - (lr(1)*hRCA(1,1)-xA(1,1)*hrg1(1))/L;
    Flfr =  m(1)*dervy_t{1}*(hRCA(1,1)*lr(1)/L/tW(1,1)+dh/tW(1,1)*crollA(1,1)/(crollA(1,1)+crg1)); % first axle lateral load transfer right side
    Flfl = - Flfr;
    Flrr =  m(1)*dervy_t{1}*(hrg1*lf(1)/L/tW(1,2)+dh/tW(1,2)*crg1/(crollA(1,1)+crg1));
    Flrl = - Flrr;
    
    % ungrouping the LLT of the axle-group of the first unit
    FzrLt{1,1} = Flfr;
    FzlLt{1,1} = Flfl;
    for j=2:nUnitAxles
        if ua(1,j)>0 % if axle exists
            FzrLt{1,j} = Flrr/(nAxUnit(1)-1);
            FzlLt{1,j} = Flrl/(nAxUnit(1)-1);
        end
    end
    
    % LLT for the units other than the first unit
    for i = 2:nUnits   
        hrgi = sum(hRCA(i,:).*ua(i,:))/(nAxUnit(i)); % group-axle roll center height
        crgi = sum(crollA(i,:).*ua(i,:)); % group-axle roll stiffness
        Flr(i) = m(i)*dervy_t{i}*h(i)/tW(i,1)/(1-m(i)*g/crgi*(h(i)-hrgi)); %Eq. 4.31 compendium 2018
        Fll(i) = - Flr(i);
    end
  
    % ungrouping the LLT of the axle-group of the units other than the
    % first unit
    for i=2:nUnits
        for j=1:nUnitAxles
            if ua(i,j)>0 % if axle exists
                FzrLt{i,j} = Flr(i)/nAxUnit(i);
                FzlLt{i,j} = Fll(i)/nAxUnit(i);
            end
        end
    end
    
end


%% Lateral Slip  

slips = cell(nUnits,size(derPW_tWCoor,2));
for i=1:nUnits
    for j= 1:nUnitAxles
        if ua(i,j)>0 % if axle exists 
            slips{i,j} = derPW_tWCoor{i,j}(2)/abs(derPW_tWCoor{i,j}(1)); % vy/abs(vx) of each wheel in wheel's coordinate system (right wheels in case of a two-track vehicle)
            if ~singleTrack
                  slips{i,j+nUnitAxles} = derPW_tWCoor{i,j+nUnitAxles}(2)/abs(derPW_tWCoor{i,j+nUnitAxles}(1)); % vy/abs(vx) of each left wheels in wheel's coordinate system
            end
        end
    end 
end

%% Nonlinear tyre model (a version of Pacejka model). 

Fyw = cell(nUnits,size(derPW_tWCoor,2));
if ~linearTyre
    u2 = 0.8;
    CC = 2*(1+asin(u2)/pi);
    uy0 = 0.8;
    ccyg = -0.1;
    Fywls = cell(nUnits,nUnitAxles);
    Fywrs = cell(nUnits,nUnitAxles);
    for i=1:nUnits
        for j= 1:nUnitAxles
            if ua(i,j)>0 % if axle exists 
                if sa(i,j)>0 %|| (i==3 && j==2) % if axle steerable (steering tyres have different porperties) (for dolly, second axle could also have steerable tyres properties)
                    uyg = uyg_steerWh;
                    CCy0 = CCy0_steerWh;
                else
                    uyg = uyg_normalWh;
                    CCy0 = CCy0_normalWh;
                end
                slipw = atan(slips{i,j});
                Fz0 = Fz0wNom(i,j);
                Fzw = Fz0A(i,j)/nW(i,j) * cos(pitchAng(i)) * cos(bankAng(i)); % wheel static vertical force normal to the road surface
                Fzlw = Fzw + FzlLt{i,j}/(nW(i,j)/2); % one left wheel vertical force with load transfer
                uy = uy0*(1+uyg*(Fzlw-Fz0)/Fz0);     % uyg: a variable to tune (-0.3<uyg<-0.1)
                CCy = CCy0*(1+ccyg*(Fzlw-Fz0)/Fz0);  % CCy0 variable to tune
                Fywls{i,j} = Fzlw*uy*sin(CC*atan(-CCy/CC/uy*slipw))*(nW(i,j)/2); % left wheels lateral force

                if ~singleTrack 
                    slipw = atan(slips{i,j+nUnitAxles});
                end
                Fzrw = Fzw + FzrLt{i,j}/(nW(i,j)/2); % one right wheel vertical force with load transfer
                uy = uy0*(1+uyg*(Fzrw-Fz0)/Fz0);     % uyg: a tuning variable (-0.3<uyg<-0.1)
                CCy = CCy0*(1+ccyg*(Fzrw-Fz0)/Fz0);  % CCy0 tuning variable
                Fywrs{i,j} = Fzrw*uy*sin(CC*atan(-CCy/CC/uy*slipw))*(nW(i,j)/2); % right wheels lateral force

                if singleTrack
                    Fz = (Fzlw+Fzrw)*nW(i,j)/2;
                    if includeCombinedSlip
                        Fyw{i,j} = real(sqrt(1-(Fxw(i,j)/mu(i,j)/Fz/e)^2))*(Fywls{i,j}+Fywrs{i,j}); % considering combined slip with Pythagoras model (friction circle or ellipse) in wheel frame
                    else   
                        Fyw{i,j} = (Fywls{i,j}+Fywrs{i,j}); % considering combined slip with Pythagoras model (friction circle or ellipse) in wheel frame
                    end
                else
                    Fzl = Fzlw*nW(i,j)/2;
                    Fzr = Fzrw*nW(i,j)/2;
                    if includeCombinedSlip
                        Fyw{i,j}            = real(sqrt(1-(Fxw(i,j)/mu(i,j)/Fzr/e)^2)) * Fywrs{i,j}; % considering combined slip with Pythagoras model (friction circle or ellipse) in wheel frame, right wheels
                        Fyw{i,j+nUnitAxles} = real(sqrt(1-(Fxw(i,j+nUnitAxles)/mu(i,j+nUnitAxles)/Fzl/e)^2)) * Fywls{i,j}; % considering combined slip with Pythagoras model (friction circle or ellipse) in wheel frame, left wheels
                    else
                        Fyw{i,j}            = Fywrs{i,j}; 
                        Fyw{i,j+nUnitAxles} = Fywls{i,j}; 
                    end
                end
            end
        end 
    end
end
        
%% Linear tyre model 

if linearTyre
    for i=1:nUnits
        for j= 1:nUnitAxles
            if ua(i,j)>0 % if axle exists 
                Fzw = Fz0A(i,j)/nW(i,j)* cos(pitchAng(i)) * cos(bankAng(i)); % wheel static vertical force normal to the road surface
                Fzlw = Fzw + FzlLt{i,j}/(nW(i,j)/2); % one left wheel vertical force with load transfer
                Fzrw = Fzw + FzrLt{i,j}/(nW(i,j)/2); % one right wheel vertical force with load transfer
                if singleTrack
                    Fz = (Fzlw+Fzrw)*nW(i,j)/2;
                    if includeCombinedSlip
                        Fyw{i,j} = real(sqrt(1-(Fxw(i,j)/mu(i,j)/Fz/e)^2))* (-C(i,j))*slips{i,j} * nW(i,j); % considering combined slip with Pythagoras model (friction circle or ellipse) in wheel frame. (Note that Fxw(i,j) is a side force not a single wheel force.)
                    else
                        Fyw{i,j} = -C(i,j)*slips{i,j} * nW(i,j); % considering combined slip with Pythagoras model (friction circle or ellipse) in wheel frame. (Note that Fxw(i,j) is a side force not a single wheel force.)
                    end
                else
                    Fzl = Fzlw*nW(i,j)/2; % left-side vertical force
                    Fzr = Fzrw*nW(i,j)/2; % right-side vertical force
                    if includeCombinedSlip
                        Fyw{i,j}            = real(sqrt(1-(Fxw(i,j)/mu(i,j)/Fzr/e)^2)) * (-C(i,j)) * slips{i,j} * nW(i,j)/2; % considering combined slip with Pythagoras model (friction circle or ellipse) in wheel frame, right wheels (Note that Fxw(i,j) is a side force not a single wheel force.)
                        Fyw{i,j+nUnitAxles} = real(sqrt(1-(Fxw(i,j+nUnitAxles)/mu(i,j+nUnitAxles)/Fzl/e)^2)) * (-C(i,j+nUnitAxles)) * slips{i,j+nUnitAxles} * nW(i,j)/2; % considering combined slip with Pythagoras model (friction circle or ellipse) in wheel frame, left wheels
                    else
                        Fyw{i,j}            = -C(i,j) * slips{i,j} * nW(i,j)/2; % considering combined slip with Pythagoras model (friction circle or ellipse) in wheel frame, right wheels (Note that Fxw(i,j) is a side force not a single wheel force.)
                        Fyw{i,j+nUnitAxles} = -C(i,j+nUnitAxles) * slips{i,j+nUnitAxles} * nW(i,j)/2; % considering combined slip with Pythagoras model (friction circle or ellipse) in wheel frame, left wheels
                    end
                end
            end
        end
    end
end

%% Lateral forces in global frame

F = cell(nUnits, size(derPW_tWCoor,2));
for i=1:nUnits
    for j= 1:nUnitAxles
        if ua(i,j)>0 % if axle exists 
            if sa(i,j)>0  % if axle steerable
                wheelAngle = phi{i} + del(i,j); 
            else
                wheelAngle = phi{i};
            end
            F{i,j} = [cos(wheelAngle) , -sin(wheelAngle);
                      sin(wheelAngle) , cos(wheelAngle)] * [Fxw(i,j);Fyw{i,j}];
            F{i,j} = F{i,j}(t);
            if ~singleTrack
                if sa(i,j)>0  % if axle steerable
                    wheelAngle = phi{i} + del(i,j+nUnitAxles); 
                else
                    wheelAngle = phi{i};
                end
                F{i,j+nUnitAxles} = [cos(wheelAngle) , -sin(wheelAngle);
                                     sin(wheelAngle) , cos(wheelAngle)] * [Fxw(i,j+nUnitAxles);Fyw{i,j+nUnitAxles}];
                F{i,j+nUnitAxles} = F{i,j+nUnitAxles}(t);
            end
        end
    end
end

%% Body forces: air, grade, and rolling resistance

% Unit body force
Fg = cell(nUnits,1);
for i=1:nUnits
    % unit local forces  
    if includeRoadGradeForce
        Fgx = m(i)*g*sin(pitchAng(i)); % angle positive downhill
        Fgy = m(i)*g*sin(bankAng(i));  % road banking angle around local x-axis, positive downhill at the left (driver) side
        % unit global forces
        Fg{i}= [cos(phi{i}) , -sin(phi{i});
                sin(phi{i}) , cos(phi{i})] * [Fgx;Fgy];
    else
        Fg{i}= [0;0];
    end
end

% Air resistance force
if includeAirResistance 
    Fair = [cos(phi{1}) ,  -sin(phi{1});
           sin(phi{1}) , cos(phi{1})] * [-0.5*AfcdRhoa*vx^2*sign(vx);0]; % no side air force
else
    Fair = [0;0];
end

% Rolling resistance forces 
FrollG = cell(nUnits, size(derPW_tWCoor,2));
for i=1:nUnits
    for j=1:nUnitAxles
        if ua(i,j)>0 % if axle exists 
            if sa(i,j)>0  % if axle steerable
                wheelAngle = phi{i} + del(i,j); 
            else
                wheelAngle = phi{i};
            end 
            Fzw = Fz0A(i,j)/2 * cos(pitchAng(i)) * cos(bankAng(i)); % left or right static vertical force normal to the road surface
            Fzlw = Fzw + FzlLt{i,j}; % left side vertical force with load transfer
            Fzrw = Fzw + FzrLt{i,j}; % right side vertical force with load transfer
            if singleTrack
                if includeRollingResistance
                    Frollw = [-(Fzlw+Fzrw) * fr * sign(vx);0]; % rolling resistance forces in wheel coord.
                else
                    Frollw = [0;0];
                end
                FrollG{i,j} = [cos(wheelAngle) , -sin(wheelAngle); % rolling resistance forces in global coord.
                               sin(wheelAngle) , cos(wheelAngle)] * Frollw;
            else 
                if sa(i,j)>0  % if axle steerable
                    wheelAngle = phi{i} + del(i,j); 
                else
                    wheelAngle = phi{i};
                end
                if includeRollingResistance
                    Frollrw = [-Fzrw * fr * sign(vx);0];   % rolling resistance forces in right wheel coord.
                else
                    Frollrw = [0;0];
                end
                FrollG{i,j} = [cos(wheelAngle) ,  -sin(wheelAngle);      % in global coord.
                               sin(wheelAngle) , cos(wheelAngle)] * Frollrw;

                if sa(i,j)>0  % if axle steerable
                    wheelAngle = phi{i} + del(i,j+nUnitAxles); 
                else
                    wheelAngle = phi{i};
                end
                if includeRollingResistance
                    Frolllw = [-Fzlw * fr * sign(vx);0];   % rolling resistance forces in left wheel coord.
                else
                    Frolllw = [0;0];
                end

                FrollG{i,j+nUnitAxles}  = [cos(wheelAngle) , -sin(wheelAngle); % in global coord.
                                           sin(wheelAngle) , cos(wheelAngle)] * Frolllw;
            end
        end 
    end
end
            


%% dT/dq, dT/dqdot

derT_phi1 = A*vy-B*vx+functionalDerivative(T,phi{1}); % the last term is always zero
derT_derPhi1_t = functionalDerivative(T,derPhi_t{1});
Eq_derT_phi1='derT_phi1 = A*vy-B*vx;'; % needed for printing the output

if nUnits>1
    for i=1:nUnits-1 
        derT_thet{i} = functionalDerivative(T,thet(i));
        derT_derThet_t{i} = functionalDerivative(T,derThet_t(i));
        Eq_derT_thet{i} = [['derT_thet',num2str(i),' = '],char(derT_thet{i}),';']; Eq_derT_thet{i} = erase(Eq_derT_thet{i},"(t)");
    end
end



%% d(dT/dqdot)/dt
 
derderT_derPhi1_t  = diff(derT_derPhi1_t,t);
derderT_derPhi1_t = subs(derderT_derPhi1_t,diff(vx(t),t),derVx_t(t));
derderT_derPhi1_t = subs(derderT_derPhi1_t,diff(vy(t),t),derVy_t(t));
derderT_derPhi1_t = subs(derderT_derPhi1_t,diff(phi1(t),t),derPhi_t1(t));
derderT_derPhi1_t = subs(derderT_derPhi1_t,diff(derPhi_t1(t),t),derderPhi_t1(t));
if nUnits>1
    for k = 1:nUnits-1
        derderT_derPhi1_t = subs(derderT_derPhi1_t,diff(thet(k),t),derThet_t(k));
        derderT_derPhi1_t = subs(derderT_derPhi1_t,diff(derThet_t(k),t),derderThet_t(k));
    end
end

Eq_derderT_derPhi1_t=['derderT_derPhi1_t = ',char(derderT_derPhi1_t),';'];    Eq_derderT_derPhi1_t = erase(Eq_derderT_derPhi1_t,"(t)"); % converting to string and removing '(t)'

derderT_derThet_t = cell(nUnits-1,1);
Eq_derderT_derThet_t = cell(nUnits-1,1);
if nUnits>1
    for i=1:nUnits-1 
        derderT_derThet_t{i} = diff(derT_derThet_t{i},t);
        derderT_derThet_t{i} = subs(derderT_derThet_t{i},diff(vx(t),t),derVx_t(t));
        derderT_derThet_t{i} = subs(derderT_derThet_t{i},diff(vy(t),t),derVy_t(t));
        derderT_derThet_t{i} = subs(derderT_derThet_t{i},diff(phi1(t),t),derPhi_t1(t));
        derderT_derThet_t{i} = subs(derderT_derThet_t{i},diff(derPhi_t1(t),t),derderPhi_t1(t));
        if nUnits>1
            for k = 1:nUnits-1
                derderT_derThet_t{i} = subs(derderT_derThet_t{i},diff(thet(k),t),derThet_t(k));
                derderT_derThet_t{i} = subs(derderT_derThet_t{i},diff(derThet_t(k),t),derderThet_t(k));
            end
        end
        Eq_derderT_derThet_t{i}=[['derderT_derThet_t',num2str(i),' = '],char(derderT_derThet_t{i}),';']; Eq_derderT_derThet_t{i} = erase(Eq_derderT_derThet_t{i},"(t)"); % converting to string and removing '(t)'
    end
end

%% Printing vehicle and equation information into a txt file

fid = fopen('equationsOfMotion.m', 'wt');

if singleTrack
    nTrack = 'single-track';
else
    nTrack = 'two-track';
end
if linearTyre
    tyreModel = 'linear';
else
    tyreModel = 'nonlinear model (a version of the Pocejka magic tyre model)';
end
if includeCombinedSlip
    combinedSlip = 'considering the friction ellipse combined slip model.';
else
    combinedSlip = 'considering no combined slip model.';
end
if includeLateralLoadTransfer
LLT = ' a steady state lateral load transfer assuming roll-free couplings, ';
else 
LLT = ' no lateral load transfer, ';
end
  
if includeRoadGradeForce
RG = ' body force parallel to the road surface, ';
else 
RG = ' no body force parallel to the road surface, ';
end

if includeAirResistance
AR = ' air resistance force, ';
else 
AR = ' no air resistance force, ';
end

if includeRollingResistance
RR = ' rolling resistance force.';
else 
RR = ' no rolling resistance force.';
end


thetStates = [];
derthetStates = [];
derderthetStates = [];
EqsInSpaceDomain = 0;
if ~EqsInSpaceDomain
    for i=1:nUnits-1
        thetStates = [thetStates, ', thet',num2str(i)];
    end
    for i=1:nUnits-1
        derthetStates = [derthetStates, ', derThet_t',num2str(i)];
    end
    for i=1:nUnits-1
        derderthetStates = [derderthetStates, ', derderThet_t',num2str(i)];
    end
end


% preparing for printing inputs and parameters:
inpsFx =[];
inpsDel = [];
ms = []; 
Js = [];
xAs = []; tWs = []; xCs = []; hs = []; nWs = []; hRCAs = []; crollAs = []; Fz0As = []; Fz0wNoms = [];
FzlLts = []; FzrLts = []; pitchAngs = []; bankAngs = []; mus = []; Cs = [];
for i=1:nUnits
    ms = [ms, char(m(i)), ', '];
    Js = [Js, char(J(i)), ', '];
    hs = [hs, char(h(i)), ', '];
    pitchAngs = [pitchAngs, char(pitchAng(i)), ', '];
    bankAngs = [bankAngs, char(bankAng(i)), ', '];
    xCs = [xCs,char(xC(i,1)),', '];
    xCs = [xCs,char(xC(i,2)),', '];
    for j= 1:nUnitAxles
        if ua(i,j)>0 % if axle exists 
            inpsFx = [inpsFx,char(Fxw(i,j)),', '];
            if sa(i,j)>0
                inpsDel = [inpsDel, char(del(i,j)),', '];
            end
            mus = [mus,char(mu(i,j)),', '];
            Cs = [Cs,char(C(i,j)),', '];
            
            xAs = [xAs,char(xA(i,j)),', '];
            tWs = [tWs,char(tW(i,j)),', '];
            nWs = [nWs,char(nW(i,j)),', '];
            hRCAs = [hRCAs,char(hRCA(i,j)),', '];
            crollAs = [crollAs,char(crollA(i,j)),', '];
            Fz0As = [Fz0As,char(Fz0A(i,j)),', '];
            Fz0wNoms = [Fz0wNoms,char(Fz0wNom(i,j)),', '];
%             FzlLts = [FzlLts,char(FzlLt{i,j}),', '];
%             FzrLts = [FzrLts,char(FzrLt{i,j}),', '];
        
            
            if ~singleTrack
                inpsFx = [inpsFx,char(Fxw(i,j+nUnitAxles)), ', '];
                mus = [mus,char(mu(i,j+nUnitAxles)),', '];
                Cs = [Cs,char(C(i,j+nUnitAxles)),', '];
                if sa(i,j)>0
                    inpsDel = [inpsDel,char(del(i,j+nUnitAxles)), ', '];
                end
            end
        end
    end
    end
inps = [inpsFx,inpsDel];
inps = inps(1:end-2);
inps= erase(inps,"(t)"); 

fprintf(fid, ['function res = equationsOfMotion(t,y,yp,param)','\n','\n']);
fprintf(fid, ['%% Nonlinear implicit equations of motion of a ',nTrack,', ',num2str(nUnits),' units articulated vehicle, ','\n']);
fprintf(fid, ['%% with ',num2str(sum(sum(ua))),' axles and ', num2str(sum(sum(sa))),' steerable axles.','\n']);

fprintf(fid, ['%% The tyre model is ',tyreModel,', and ', combinedSlip,'\n']);
fprintf(fid, ['%% The model includes: \n']);
fprintf(fid, ['%%',LLT,'\n']);
fprintf(fid, ['%%',RG,'\n']);
fprintf(fid, ['%%',AR,'\n']);
fprintf(fid, ['%%',RR,'\n','\n']);
fprintf(fid, ['%% States are x = [vx, vy, phi1, derPhi_t1',thetStates,derthetStates,'].','\n']);
fprintf(fid, ['%% States derivatives are x = [derVx_t, derVy_t, derPhi_t1',derthetStates,derderthetStates,'].','\n']);
fprintf(fid, ['%% Inputs are u = [',inps,'].','\n']);
if singleTrack
    fprintf(fid, ['%% where, Fxwi_j is the longitudinal axle force in axle wheel frame, i: unit index, j: axle index;','\n']);
    fprintf(fid, ['%% and, similarly, deli_j denote steering angles. ','\n','\n']);
else
    fprintf(fid, ['%% where, Fxwi_j is the longitudinal wheel force in wheel frame, i: unit index, j: wheel index (left-side wheels are indexed after right-side wheels);','\n']);
    fprintf(fid, ['%% and, similarly, deli_j denote steering angles. ','\n','\n']);
end

fprintf(fid, ['%% User-defined parameters: ','\n','\n']);
fprintf(fid, ['%% nUnits, %% Number of vehicle units ','\n']);
fprintf(fid, ['%% ', ms, '                     %% Unit masses ','\n']);
fprintf(fid, ['%% ', Js, '                     %% Unit inertias ','\n']);
fprintf(fid, ['%% ', xAs, '               %% Local x position of the axle w.r.t. unit COG ','\n']);
fprintf(fid, ['%% ', xCs, '               %% Unit coupling local x position w.r.t. unit COG','\n']);
if ~singleTrack || includeLateralLoadTransfer
    fprintf(fid, ['%% ', tWs, '               %% Axles track width ','\n']);
end
fprintf(fid, ['%% ', nWs, '               %% Number of wheels per axle ','\n']);
if includeLateralLoadTransfer
    fprintf(fid, ['%% ', hs, '                     %% Unit COG height ','\n']);
    fprintf(fid, ['%% ', hRCAs, '             %% Axle roll center height','\n']);
    fprintf(fid, ['%% ', crollAs, '           %% Axle rolling stiffness ','\n']);
end

if ~linearTyre
    fprintf(fid, ['%% ', Fz0wNoms, '             %% Tyre nominal vertical load ','\n']);
end
if includeCombinedSlip || ~linearTyre || includeRollingResistance
    fprintf(fid, ['%% ', Fz0As, '             %% Axle static vertical load ','\n']);
end
if includeRoadGradeForce || includeCombinedSlip
    fprintf(fid, ['%% ', pitchAngs, '          %% Road pitch angle at cog of each unit (positive downhill)','\n']);
    fprintf(fid, ['%% ', bankAngs, '           %% Road banking angle around local x-axis, positive downhill at the left (driver) side ','\n']);
end
if includeCombinedSlip
    fprintf(fid, ['%% e,                      %% friction circle scale factor (if e is not 1 then the circle changes to an ellipse)','\n']);
end
fprintf(fid, ['%% g, AfcdRhoa, fr,           %% Gravitational constant, air resistance coefficients, rolling resistance coefficient (if included in the model)','\n']);
if includeCombinedSlip || ~linearTyre
    fprintf(fid, ['%% ', mus, '               %% Road friction coefficient (can be different for different axles and sides depending on road position)','\n']);
end
if linearTyre 
    fprintf(fid, ['%% ', Cs, '            %%  A single wheel cornering stiffness (can be different for different axles and sides)','\n','\n']);
else
    fprintf(fid, ['%% uyg_steerWh, CCy0_steerWh, uyg_normalWh, CCy0_normalWh, %% Nonlinear tyre model tuning parameters, different for steerable wheels and unsteerable (normal) wheels','\n','\n']);
end


fprintf(fid, ['%% States: ','\n','\n']);
fprintf(fid, ['phi1 = y(3);','\n']);
for i = 1:nUnits-1
    fprintf(fid, ['thet',num2str(i),' = y(',num2str(i+3),');','\n']);
end
fprintf(fid, ['vx = y(',num2str(nUnits+3),');','\n']);
fprintf(fid, ['vy = y(',num2str(nUnits+4),');','\n']);
fprintf(fid, ['derPhi_t1 = y(',num2str(nUnits+5),');','\n']);
for i = 1:nUnits-1
    fprintf(fid, ['derThet_t',num2str(i),' = y(',num2str(i+nUnits+5),');','\n']);
end
fprintf(fid, '\n');

fprintf(fid, ['%% State derivatives: ','\n','\n']);
fprintf(fid, ['derVx_t = yp(',num2str(nUnits+3),');','\n']);
fprintf(fid, ['derVy_t = yp(',num2str(nUnits+4),');','\n']);
fprintf(fid, ['derderPhi_t1 = yp(',num2str(nUnits+5),');','\n']);
for i = 1:nUnits-1
    fprintf(fid, ['derderThet_t',num2str(i),' = yp(',num2str(i+nUnits+5),');','\n']);
end
fprintf(fid, '\n');

fprintf(fid, ['%% Inputs: ','\n','\n']);

for i=1:nUnits
    for j= 1:nUnitAxles
        if ua(i,j)>0 % if axle exists 
            fprintf(fid, [erase(char(Fxw(i,j)),"(t)"), '= param.Fxw(',num2str(i),',',num2str(j),');  %% Input force','\n']);
            if sa(i,j)>0
                fprintf(fid, [erase(char(del(i,j)),"(t)"), '= param.del(',num2str(i),',',num2str(j),');  %% Input steering angle','\n']);
            end

            if ~singleTrack
                jj = j+nUnitAxles;
                fprintf(fid, [erase(char(Fxw(i,jj)),"(t)"), '= param.Fxw(',num2str(i),',',num2str(jj),');  %% Input force, left wheel','\n']);
                if sa(i,j)>0
                    fprintf(fid, [erase(char(del(i,jj)),"(t)"), '= param.del(',num2str(i),',',num2str(jj),');  %% Input steering angle, left wheel','\n']);
                end
            end
        end
    end 
    fprintf(fid,'\n','\n');
end

fprintf(fid,'%% A very simple proportional speed controller (controlling the second axle longitudinal force):\n');

if singleTrack 
    fprintf(fid,'Fxw1_2=(param.vxRef(floor(t*param.sampleRate)+1)-vx)*(sum(param.m));\n \n');
else
    fprintf(fid,'Fxw1_2=(param.vxRef(floor(t*param.sampleRate)+1)-vx)*0.5*(sum(param.m));\n');
    fprintf(fid,['Fxw1_',num2str(2+nUnitAxles),'=(param.vxRef(floor(t*param.sampleRate)+1)-vx)*0.5*(sum(param.m));\n \n']);
end


fprintf(fid,'%% Reading the steering input of the first axle: \n');

if singleTrack 
    fprintf(fid,'if t <1/param.sampleRate \n');
    fprintf(fid,'    del1_1 =0;\n');
    fprintf(fid,'else \n');
    fprintf(fid,'    del1_1 = param.del1Ref(floor(t*param.sampleRate));\n');
    fprintf(fid,'end \n\n');
else
    fprintf(fid,'if t <1/param.sampleRate \n');
    fprintf(fid,'    del1_1 =0;\n');
    fprintf(fid,['    del1_',num2str(1+nUnitAxles),' =0;\n']);
    fprintf(fid,'else \n');
    fprintf(fid,'    del1_1 = param.del1_1(floor(t*param.sampleRate));\n');
    fprintf(fid,['    del1_',num2str(1+nUnitAxles),' = param.del1_1(floor(t*param.sampleRate));\n']);
    fprintf(fid,'end \n\n');
end



fprintf(fid, ['%% Reading param to set the user-defined parameters: ','\n','\n']);
fprintf(fid, ['nUnits   =',num2str(nUnits),';           %% Number of vehicle units','\n']);
for i=1:nUnits
    fprintf(fid,['%% Unit ',num2str(i),'\n']);
    fprintf(fid, [char(m(i)),'       = param.m(',num2str(i),');             %% Unit mass (kg)','\n']);
    fprintf(fid, [char(J(i)),'       = param.J(',num2str(i),');             %% Unit inertia (km m^2)','\n']);
    if includeLateralLoadTransfer
        fprintf(fid, [char(h(i)),'       = param.h(',num2str(i),');             %% Unit COG height (m)','\n']);
    end
    if includeRoadGradeForce || includeCombinedSlip
        fprintf(fid, [char(pitchAng(i)),'= param.pitchAng(',num2str(i),');      %% Road pitch angle at cog of each unit (positive downhill) (rad)','\n']);
        fprintf(fid, [char(bankAng(i)),' = param.bankAng(',num2str(i),');       %% Road banking angle around local x-axis, positive downhill at the left side (rad)','\n']);
    end
    fprintf(fid, [char(xC(i,1)),'      = param.xC(',num2str(i),',',num2str(1),');          %% Unit coupling local x position w.r.t. unit COG','\n']);
    fprintf(fid, [char(xC(i,2)),'      = param.xC(',num2str(i),',',num2str(2),');          %% Unit coupling local x position w.r.t. unit COG','\n']);

    for j= 1:nUnitAxles
        if ua(i,j)>0 % if axle exists 
            if includeCombinedSlip || ~linearTyre
                fprintf(fid, [char(mu(i,j)),'      = param.mu(',num2str(i),',',num2str(j),');          %% Road friction','\n']);
            end
            if linearTyre
                fprintf(fid, [char(C(i,j)),'       = param.C(',num2str(i),',',num2str(j),');           %% Cornering stiffness (N)','\n']);
            end
            fprintf(fid, [char(xA(i,j)),'      = param.xA(',num2str(i),',',num2str(j),');          %% Local x position of the axle w.r.t. unit COG (m)','\n']);
            fprintf(fid, [char(nW(i,j)),'      = param.nW(',num2str(i),',',num2str(j),');          %% Number of wheels per axle','\n']);
            if ~singleTrack || includeLateralLoadTransfer
                fprintf(fid, [char(tW(i,j)),'      = param.tW(',num2str(i),',',num2str(j),');          %% Axles track width (m)','\n']);
            end
            if includeLateralLoadTransfer
                fprintf(fid, [char(hRCA(i,j)),'    = param.hRCA(',num2str(i),',',num2str(j),');        %% Axle roll center height (m)','\n']);
                fprintf(fid, [char(crollA(i,j)),'  = param.crollA(',num2str(i),',',num2str(j),');      %% Axle rolling stiffness (N m/rad) ','\n']);
            end
            if ~linearTyre 
                fprintf(fid, [char(Fz0wNom(i,j)),'    = param.Fz0wNom(',num2str(i),',',num2str(j),');        %% Wheel nominal vertical load (N)','\n']);
            end
            if includeCombinedSlip || ~linearTyre || includeRollingResistance
                fprintf(fid, [char(Fz0A(i,j)),'    = param.Fz0A(',num2str(i),',',num2str(j),');        %% Axle static vertical load (N)','\n']);
            end
            if ~singleTrack
                jj = j+nUnitAxles;
                if includeCombinedSlip || ~linearTyre
                    fprintf(fid,     [char(mu(i,jj)),'      = param.mu(',num2str(i),',',num2str(jj),');         %% Road friction, left wheel','\n']);
                end
                if linearTyre
                    fprintf(fid, [char(C(i,jj)),'   = param.C(',num2str(i),',',num2str(jj),');              %% Cornering stiffness, left wheel (N)','\n']);
                end
            end
        end
    end 
    fprintf(fid,'\n');
end
if includeCombinedSlip
    fprintf(fid, ['e     = param.e;             %% friction ellipse scale factor (if e is 1 then the ellipse becomes a circle)','\n']);
end
if includeRoadGradeForce || includeLateralLoadTransfer
    fprintf(fid, ['g        = param.g;                %% Gravitational acceleration (m/s^2)','\n']);
end
if includeAirResistance
    fprintf(fid, ['AfcdRhoa = param.AfcdRhoa;         %% air resistance coefficients (N/(m/s)^2)','\n']);
end
if includeRollingResistance
    fprintf(fid, ['fr       = param.fr;               %% rolling resistance coefficient','\n']);
end

if ~linearTyre 
    fprintf(fid, ['uyg_steerWh   = param.uyg_steerWh; %% Nonlinear tyre model tuning parameters, different for steerable wheels and unsteerable (normal) wheels','\n']);
    fprintf(fid, ['CCy0_steerWh  = param.CCy0_steerWh;','\n']);
    fprintf(fid, ['uyg_normalWh  = param.uyg_normalWh;','\n']);
    fprintf(fid, ['CCy0_normalWh = param.CCy0_normalWh;','\n']);
end
fprintf(fid,'\n');
fprintf(fid,'%%%% Equations: \n');
%% Extracting the force elements from equations so that they will be evaluated only once during one time step integration

if extractForceElementsFromEqs
    
    % The below assignments are needed to make exctraction of the vector
    % elements possible, otherwise F1f(1) gives F1f(t=1) .
    for i=1:nUnits
        if includeRoadGradeForce
            Fg{i,1} = Fg{i,1}(t);
        else
            Fg{i,1} = ['0';'0'];
        end
        for j=1:nUnitAxles
            if ua(i,j)>0 % if axle exists 
                if includeRollingResistance
                    FrollG{i,j} = FrollG{i,j}(t);
                else
                    FrollG{i,j} = ['0';'0'];
                end
                if ~singleTrack
                    if includeRollingResistance
                        FrollG{i,j+nUnitAxles} = FrollG{i,j+nUnitAxles}(t);
                    else
                        FrollG{i,j+nUnitAxles} =['0';'0'];
                    end
                end
            end
        end
    end
    if includeAirResistance
        Fair = Fair(t);
    else
        Fair = ['0';'0'];
    end
    
    
    % Creating strings of force elements experessions needed for printing to a txt file
    
    fe_Fairx = ['Fair1 = ',char(Fair(1)),';']; fe_Fairx = erase(fe_Fairx,"(t)"); % converting to string and removing '(t)'
    fe_Fairy = ['Fair2 = ',char(Fair(2)),';']; fe_Fairy = erase(fe_Fairy,"(t)"); % converting to string and removing '(t)'
    fprintf(fid, [fe_Fairx,'\n', fe_Fairy,'\n']);

    for i=1:nUnits
        fe_Fgx{i,1} = ['Fg',num2str(i),'1 = ',char(Fg{i,1}(1)),';']; fe_Fgx{i,1} = erase(fe_Fgx{i,1},"(t)"); % converting to string and removing '(t)'
        fe_Fgy{i,1} = ['Fg',num2str(i),'2 = ',char(Fg{i,1}(2)),';']; fe_Fgy{i,1} = erase(fe_Fgy{i,1},"(t)"); % converting to string and removing '(t)'
        fprintf(fid, [fe_Fgx{i,1},'\n', fe_Fgy{i,1},'\n']);
        for j=1:nUnitAxles
            if ua(i,j)>0 % if axle exists 
                fe_Fx{i,j} = ['F',num2str(i),num2str(j),'1 = ',char(F{i,j}(1)),';']; fe_Fx{i,j} = erase(fe_Fx{i,j},"(t)"); % converting to string and removing '(t)'
                fe_Fy{i,j} = ['F',num2str(i),num2str(j),'2 = ',char(F{i,j}(2)),';']; fe_Fy{i,j} = erase(fe_Fy{i,j},"(t)"); % converting to string and removing '(t)'

                fe_FrollGx{i,j} = ['FrollG',num2str(i),num2str(j),'1 = ',char(FrollG{i,j}(1)),';']; fe_FrollGx{i,j} = erase(fe_FrollGx{i,j},"(t)"); % converting to string and removing '(t)'
                fe_FrollGy{i,j} = ['FrollG',num2str(i),num2str(j),'2 = ',char(FrollG{i,j}(2)),';']; fe_FrollGy{i,j} = erase(fe_FrollGy{i,j},"(t)"); % converting to string and removing '(t)'
                
                fprintf(fid, [fe_Fx{i,j},'\n', fe_Fy{i,j},'\n',fe_FrollGx{i,j},'\n', fe_FrollGy{i,j},'\n']);
                
                if ~singleTrack
                    jj = j + nUnitAxles;
                    fe_Fx{i,jj} = ['F',num2str(i),num2str(jj),'1 = ',char(F{i,jj}(1)),';']; fe_Fx{i,jj} = erase(fe_Fx{i,jj},"(t)"); % converting to string and removing '(t)'
                    fe_Fy{i,jj} = ['F',num2str(i),num2str(jj),'2 = ',char(F{i,jj}(2)),';']; fe_Fy{i,jj} = erase(fe_Fy{i,jj},"(t)"); % converting to string and removing '(t)'

                    fe_FrollGx{i,jj} = ['FrollG',num2str(i),num2str(jj),'1 = ',char(FrollG{i,jj}(1)),';']; fe_FrollGx{i,jj} = erase(fe_FrollGx{i,jj},"(t)"); % converting to string and removing '(t)'
                    fe_FrollGy{i,jj} = ['FrollG',num2str(i),num2str(jj),'2 = ',char(FrollG{i,jj}(2)),';']; fe_FrollGy{i,jj} = erase(fe_FrollGy{i,jj},"(t)"); % converting to string and removing '(t)'
                    
                    fprintf(fid, [fe_Fx{i,jj},'\n', fe_Fy{i,jj},'\n',fe_FrollGx{i,jj},'\n', fe_FrollGy{i,jj},'\n']);
                end
            end
        end
    end
    
    
    % since the force elements were already printed to the file, they should be
    % redefined as symbols:
    Fg = cell(nUnits, 1);
    Fs = cell(nUnits, size(F,2));
    FrollGs = cell(nUnits, size(F,2));
    for i=1:nUnits
        Fg{i,1} = sym(['Fg',num2str(i)],[2  1]); 
        for j=1:nUnitAxles
            if ua(i,j)>0 % if axle exists 
                F{i,j} = sym(['F',num2str(i),num2str(j)],[2  1]); 
                FrollG{i,j} = sym(['FrollG',num2str(i),num2str(j)],[2  1]);
                if ~singleTrack
                    F{i,j+nUnitAxles} = sym(['F',num2str(i),num2str(j+nUnitAxles)],[2  1]); 
                    FrollG{i,j+nUnitAxles} = sym(['FrollG',num2str(i),num2str(j+nUnitAxles)],[2  1]);
                end
            end
        end
    end
    Fair = sym('Fair',[2 1]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generalized forces: lateral add cog forces

% The below assignments are needed to make extraction of the vector
% elements possible, otherwise p1wf(1) gives p1wf(t=1) .
for i=1:nUnits
    pcog{i,1} = pcog{i,1}(t);
    for j=1:nUnitAxles
        if ua(i,j)>0 % if axle exists 
            pW{i,j} = pW{i,j}(t);
            if ~singleTrack
                pW{i,j+nUnitAxles} = pW{i,j+nUnitAxles}(t);
            end
        end
    end
end
pair = pair(t);

% functionalDerivative used instead of diff, because diff cannot handle
% derivative of a function with respect to another function e.g. X(t).
% functionalDerivative, however, cannot handle vector functions. 
QX = Fair.'*[functionalDerivative(pair(1),X);functionalDerivative(pair(2),X)]; % contribution of air resistance force
QY = Fair.'*[functionalDerivative(pair(1),Y);functionalDerivative(pair(2),Y)]; % contribution of air resistance force
Qphi1 = Fair.'*[functionalDerivative(pair(1),phi{1});functionalDerivative(pair(2),phi{1})]; % contribution of air resistance force
if nUnits>1
    for k=1:nUnits-1
        Qthet{k} = Fair.'*[functionalDerivative(pair(1),thet(k));functionalDerivative(pair(2),thet(k))]; % contribution of air resistance force
    end
end
for i=1:nUnits
    QX = QX + Fg{i,1}.'*[functionalDerivative(pcog{i,1}(1),X);functionalDerivative(pcog{i,1}(2),X)]; % contribution of body forces
    QY = QY + Fg{i,1}.'*[functionalDerivative(pcog{i,1}(1),Y);functionalDerivative(pcog{i,1}(2),Y)]; % contribution of body forces
    Qphi1 = Qphi1 + Fg{i,1}.'*[functionalDerivative(pcog{i,1}(1),phi{1});functionalDerivative(pcog{i,1}(2),phi{1})]; % contribution of body forces
    
    if nUnits>1
        for k=1:nUnits-1
            Qthet{k} = Qthet{k} + Fg{i,1}.'*[functionalDerivative(pcog{i,1}(1),thet(k));functionalDerivative(pcog{i,1}(2),thet(k))]; % contribution of body forces
        end
    end

    for j= 1:nUnitAxles
        if ua(i,j)>0 % if axle exists 
            QX = QX + F{i,j}.'*[functionalDerivative(pW{i,j}(1),X);functionalDerivative(pW{i,j}(2),X)]; % contribution of tyre forces
            QX = QX + FrollG{i,j}.'*[functionalDerivative(pW{i,j}(1),X);functionalDerivative(pW{i,j}(2),X)]; % contribution of rolling resistance forces
            
            QY = QY + F{i,j}.'*[functionalDerivative(pW{i,j}(1),Y);functionalDerivative(pW{i,j}(2),Y)]; % contribution of tyre forces
            QY = QY + FrollG{i,j}.'*[functionalDerivative(pW{i,j}(1),Y);functionalDerivative(pW{i,j}(2),Y)]; % contribution of rolling resistance forces

            Qphi1 = Qphi1 + F{i,j}.'*[functionalDerivative(pW{i,j}(1),phi{1});functionalDerivative(pW{i,j}(2),phi{1})]; % contribution of tyre forces
            Qphi1 = Qphi1 + FrollG{i,j}.'*[functionalDerivative(pW{i,j}(1),phi{1});functionalDerivative(pW{i,j}(2),phi{1})]; % contribution of rolling resistance forces
            
            if nUnits>1
                for k=1:nUnits-1
                    Qthet{k} = Qthet{k} + F{i,j}.'*[functionalDerivative(pW{i,j}(1),thet(k));functionalDerivative(pW{i,j}(2),thet(k))]; % contribution of tyre forces
                    Qthet{k} = Qthet{k} + FrollG{i,j}.'*[functionalDerivative(pW{i,j}(1),thet(k));functionalDerivative(pW{i,j}(2),thet(k))]; % contribution of rolling resistance forces
                end
            end

            if ~singleTrack
                QX = QX + F{i,j+nUnitAxles}.'*[functionalDerivative(pW{i,j+nUnitAxles}(1),X);functionalDerivative(pW{i,j+nUnitAxles}(2),X)]; % contribution of tyre forces (left wheels)
                QX = QX + FrollG{i,j+nUnitAxles}.'*[functionalDerivative(pW{i,j+nUnitAxles}(1),X);functionalDerivative(pW{i,j+nUnitAxles}(2),X)]; % contribution of rolling resistance forces (left wheels)
            
                QY = QY + F{i,j+nUnitAxles}.'*[functionalDerivative(pW{i,j+nUnitAxles}(1),Y);functionalDerivative(pW{i,j+nUnitAxles}(2),Y)]; % contribution of tyre forces (left wheels)
                QY = QY + FrollG{i,j+nUnitAxles}.'*[functionalDerivative(pW{i,j+nUnitAxles}(1),Y);functionalDerivative(pW{i,j+nUnitAxles}(2),Y)]; % contribution of rolling resistance forces (left wheels)

                Qphi1 = Qphi1 + F{i,j+nUnitAxles}.'*[functionalDerivative(pW{i,j+nUnitAxles}(1),phi{1});functionalDerivative(pW{i,j+nUnitAxles}(2),phi{1})]; % contribution of tyre forces (left wheels)
                Qphi1 = Qphi1 + FrollG{i,j+nUnitAxles}.'*[functionalDerivative(pW{i,j+nUnitAxles}(1),phi{1});functionalDerivative(pW{i,j+nUnitAxles}(2),phi{1})]; % contribution of rolling resistance forces (left wheels)
                
                if nUnits>1
                    for k=1:nUnits-1
                        Qthet{k} = Qthet{k} + F{i,j+nUnitAxles}.'*[functionalDerivative(pW{i,j+nUnitAxles}(1),thet(k));functionalDerivative(pW{i,j+nUnitAxles}(2),thet(k))]; % contribution of tyre forces
                        Qthet{k} = Qthet{k} + FrollG{i,j+nUnitAxles}.'*[functionalDerivative(pW{i,j+nUnitAxles}(1),thet(k));functionalDerivative(pW{i,j+nUnitAxles}(2),thet(k))]; % contribution of rolling resistance forces
                    end
                end
                
            end
        end
    end
end


% converting the QX and QY to quasi-static (first unit local coord.)
Qxy = [cos(phi{1}) , sin(phi{1});
       -sin(phi{1}) , cos(phi{1})] * [QX; QY];
Qxy = Qxy(t);
Qx = Qxy(1); Qy = Qxy(2);


Eq_Qx=['Qx = ',char(Qx),';']; Eq_Qx = erase(Eq_Qx,"(t)"); % converting to string and removing '(t)'
Eq_Qy=['Qy = ',char(Qy),';']; Eq_Qy = erase(Eq_Qy,"(t)");
Eq_Qphi1=['Qphi1 = ',char(Qphi1),';'];   Eq_Qphi1 = erase(Eq_Qphi1,"(t)");
if nUnits>1
    for i=1:nUnits-1
        Eq_Qthet{i}=[['Qthet',num2str(i),' = '],char(Qthet{i}),';']; Eq_Qthet{i} = erase(Eq_Qthet{i},"(t)");
    end
end

%% Print the equations into a txt file

fprintf(fid, [Eq_A,'\n', Eq_B,'\n',Eq_derA_t,'\n',Eq_derB_t,'\n',Eq_Qx,'\n',Eq_Qy,'\n',Eq_Qphi1,'\n',Eq_derT_phi1,'\n',Eq_derderT_derPhi1_t,'\n']);

if nUnits>1
    for i=1:nUnits-1
        fprintf(fid,[Eq_Qthet{i},'\n', Eq_derT_thet{i},'\n', Eq_derderT_derThet_t{i},'\n']); 
    end
end

% printing equation of motion:
fprintf(fid,['\n','residualVal = zeros(nUnits+2,1);','\n']);
fprintf(fid,['residualVal(1) = derA_t - B*derPhi_t1 - Qx;','\n']);
fprintf(fid,['residualVal(2) = derB_t + A*derPhi_t1 - Qy;','\n']);
fprintf(fid,['residualVal(3) = derderT_derPhi1_t - derT_phi1 - Qphi1;','\n']);
if nUnits>1
    for i=1:nUnits-1
        fprintf(fid, ['residualVal(',num2str(i+3),') = derderT_derThet_t',num2str(i),' - derT_thet',num2str(i),' - Qthet',num2str(i),';','\n']);
    end
end

fprintf(fid,['res=zeros(',num2str(2*(nUnits+2)),',1);','\n']);
fprintf(fid,['res(1:',num2str(nUnits+2),',1)=y(',num2str(nUnits+3),':',num2str(2*(nUnits+2)),')-yp(1:',num2str(nUnits+2),');','\n']);
fprintf(fid,['res(',num2str(nUnits+3),':',num2str(2*(nUnits+2)),',1) = residualVal;','\n']);

fclose(fid);

disp('Done! ')
disp('The equations are printed in file ''equationsOfMotion.m''. ')
disp('If you generated the equations for the default vehicle (the single-track 4-unit 6-axle vehicle) run the file ''vehicleSimulation.m'' to simulate the model. ')
disp('Otherwise, the parameters and inputs defined in the file ''parameters.m'' need to be updated as specified in the file ''equationsOfMotion.m''. ')
disp('Also, the solver settings and states in the file ''vehicleSimulation.m'' need to be updated based on the new model states.')


This code generates all differential equations describing the dynamic motion of a
multi-trailer vehicle, for simulation, optimization, and optimal control purposes. 
The code was generated using MATLAB R2019b, and works in this or higher versions of 
MATLAB.

Please cite the following paper when referring to this code:
T. Ghandriz, B. Jacobson, P. Nilsson, and L. Laine, “Computationally efficient
nonlinear one- and two-track models for multi-trailer road vehicles,”
IEEE Access, vol. -, pp. –, 2020.

Explanation of the vehicle model and equations can be found in the above paper. 

The first version of the code includes three Matlab files:

generateEqs.m,
parameters.m,
vehicleSimulation.m.

The file "generateEqs.m" generates the equations of motion, based on the  
user-defined vehicle, in the form of implicit ODEs, and stores them in Matlab 
function "equationsOfMotion.m", which is created automatically. The Matlab function 
"equationsOfMotion.m" is constructed in such a way that it can be directly 
called by the Matlab DAE solver, i.e., ode15i, to simulate the vehicle motion. 
In order to facilitate the simulation, we added a 
simple proportional speed controller that controls the forces of the second 
axle of the first unit to track a reference desired speed. Moreover, the 
steering input of the first axle reads the user-defined steering angle 
reference trajectory. 

As a simulation example, we provided files ``parameters.m'' and ``vehicleSimulation.m''.  
The file ``parameters.m'' defines parameters of the vehicle, and the file 
``vehicleSimulation.m'' integrates the ODEs and shows the result. These two files are 
valid for a single-track 4-unit 6-axle vehicle, and need to be updated by the user for 
simulating a different vehicle.  

In the current version of the code, a user have options to create a model including:
- any number of vehicle units;
- any number of axles in each vehicle unit;
- a single-track vehicle or a two-track vehicle;
- linear or nonlinear tyre models;
- inclusion of lateral load transfer;
- inclusion of combined slip tyre model;
- inclusion of the force caused by the road grade, and air and rolling resistance forces.
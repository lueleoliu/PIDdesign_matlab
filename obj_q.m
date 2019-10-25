function y = obj_q(x)

global QPID_type QFilter sys

Kp = x(1);
Ti = x(2);
Kd = 0;
Td = 0;

QPID = ss(tf([Kd + Kp*Td, Kp + Kp/Ti*Td, Kp/Ti], [Td, 1, 0]));
QPID = QPID*QFilter; 
QPID.inputname = 'Measured generator speed';
QPID.outputname = 'Collective generator torque demand';

% Assign Loop %

Controller = QPID;
ControllerType = QPID_type;

sysd = convertSsTime(sys,ControllerType.TimeStep,'ZOH','Plant');
sysdol = getSiso(sysd,Controller.outputname,Controller.inputname);
sysdol.outputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep);
sysdol.inputdelay = round((ControllerType.InputDelay + ControllerType.OutputDelay)/ControllerType.TimeStep)-sysdol.outputdelay;
Controller = c2d(Controller,ControllerType.TimeStep,'tustin');   

Controller.inputname = QPID.inputname;
Controller.outputname = QPID.outputname;

ControllerDelay = Controller;
ControllerDelay.inputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep); % add delays to the controller
ControllerDelay.outputdelay = round((ControllerType.InputDelay + ControllerType.OutputDelay)/ControllerType.TimeStep)-sysdol.outputdelay; % add delays to the controller

% The closed loop system %
syscl = addController(sys, ControllerDelay, ControllerType);
[y_sam,t_sam] = step(getSiso(syscl,'Collective wind speed','Generator torque'),100);
Step1 = stepinfo(y_sam,t_sam);
y(1) = Step1.SettlingTime;
y(2) = Step1.RiseTime;
    






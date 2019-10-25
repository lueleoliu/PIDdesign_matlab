function [C, Ceq] = nonl_q(x)

global QPID_type QFilter sys base_set

Kp = x(1);
Ti = x(2);
Kd = 0;
Td = 0;

QPID = ss(tf([Kd + Kp*Td, Kp + Kp/Ti*Td, Kp/Ti], [Td, 1, 0]));
QPID = QPID*QFilter; 
QPID.inputname = 'Measured generator speed';
QPID.outputname ='Collective generator torque demand';

% Assign Loop %

Controller = QPID;
ControllerType = QPID_type;

sysd = convertSsTime(sys,ControllerType.TimeStep,'ZOH','Plant');
sysdol = getSiso(sysd,Controller.outputname,Controller.inputname);
sysdol.outputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep);
sysdol.inputdelay = round((ControllerType.InputDelay + ControllerType.OutputDelay)/ControllerType.TimeStep)-sysdol.outputdelay;
Controller = c2d(Controller,ControllerType.TimeStep,'tustin');   

% The open loop system %
sysol = -Controller*sysdol;

pm_base = base_set(1,1);
wpm_base = base_set(1,2);

[~,Pm,~,Wpm] = margin(sysol);

if ~isfinite(Pm)|| isnan(Pm)
    Pm = 0;
end

if ~isfinite(Wpm)|| isnan(Wpm)
    Wpm = 0;
end

C(1) = pm_base - Pm;
C(2) = wpm_base - Wpm;

Ceq = [];


    






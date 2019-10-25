function [C, Ceq] = nonl_p(x)

global PPID_type PFilter sys base_set

Kp = x(1);
Ti = x(2);
Kd = x(3);
Td = x(4);

PPID = ss(tf([Kd + Kp*Td, Kp + Kp/Ti*Td, Kp/Ti], [Td, 1, 0]));
PPID = PPID*PFilter; 
PPID.inputname = 'Measured generator speed';
PPID.outputname ='Collective pitch angle demand';

% Assign Loop %

Controller = PPID;
ControllerType = PPID_type;

sysd = convertSsTime(sys,ControllerType.TimeStep,'ZOH','Plant');
sysdol = getSiso(sysd,Controller.outputname,Controller.inputname);
sysdol.outputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep);
sysdol.inputdelay = round((ControllerType.InputDelay + ControllerType.OutputDelay)/ControllerType.TimeStep)-sysdol.outputdelay;
Controller = c2d(Controller,ControllerType.TimeStep,'tustin');   

% The open loop system %
sysol = -Controller*sysdol;

[Gm,Pm,~,Wpm] = margin(sysol);

pm_base = base_set(1,1);
wpm_base = base_set(1,2);
gm_base = base_set(1,3);

if ~isfinite(Pm)|| isnan(Pm)
    Pm = 0;
end

if ~isfinite(Wpm)|| isnan(Wpm)
    Wpm = 0;
end

if ~isfinite(Gm)|| isnan(Gm)
    Gm = 0.1;
end


C(1) = pm_base - Pm;
C(2) = wpm_base - Wpm;
C(3) = gm_base - 20 * log10(Gm);

Ceq = [];


    






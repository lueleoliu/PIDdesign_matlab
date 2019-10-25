first_info = textread('Result.txt');

first_kp = first_info(1,1);
first_ti = first_info(1,2);
first_kd = first_info(1,3);
first_td = first_info(1,4);
first_ts = first_info(1,5);
first_tr = first_info(1,6);

load('linmod1.mat')

%% Initialise %%

% Init Filters
QFilter = 1;
PFilter = 1;
NFilter = 1;
  
% Time Steps
TsQ = 0.02;	%Torque controller timestep
TsP = 0.02;	%Pitch controller timestep
TsC = 0.02; %Controller platform timestep

% Delays
DelayW = 0.02;	%Speed input delay
DelayA = 0.02;	%Nacelle-x acceleration input delay
DelayL = 0.02;	%Load measurements input delay
DelayQ = 0.06;	%Torque output delay
DelayP = 0.08;	%Pitch output delay

% Define controller properties
CP_type.TimeStep = TsQ;
CP_type.InputDelay = DelayW;
CP_type.OutputDelay = TsC + DelayQ;
CP_type.Name = 'Constant power term';
CP_type.Prefix = 'CP';

QPID_type.TimeStep = TsQ;
QPID_type.InputDelay = DelayW;
QPID_type.OutputDelay = TsC + DelayQ;
QPID_type.Name = 'Torque-speed PID';
QPID_type.Prefix = 'QPID';

PPID_type.TimeStep = TsP;
PPID_type.InputDelay = DelayW;
PPID_type.OutputDelay = TsC + DelayP;
PPID_type.Name = 'Pitch-speed PID';
PPID_type.Prefix = 'PPID';

NAF_type.TimeStep = TsP;
NAF_type.InputDelay = DelayA;
NAF_type.OutputDelay = TsC + DelayP;
NAF_type.Name = 'Nacelle Acceleration Feedback';
NAF_type.Prefix = 'NAF';

%% Load Filters %%
filters_p = xlsread('Filters.xlsx', 'Pitch');

w = filters_p(1,3);  z= filters_p(1,4);
num = w^2;
den = [1 2*z*w w^2];
Plpf = tf(num,den);
PFilter = PFilter*ss(Plpf);

for i = 2:length(filters_p) 

    wn = filters_p(i,1); dn = filters_p(i,2); wd = filters_p(i,3); dd = filters_p(i,4);
    num = [1/wn.^2 2*dn/wn 1];
    den = [1/wd.^2 2*dd/wd 1];
    n = ss(tf(num,den));
    PFilter = PFilter*n;

end

filters_q = xlsread('Filters.xlsx', 'Torque');%%%%%%%%%%%%%%%%

w = filters_q(1,3);  z= filters_q(1,4);
num = w^2;
den = [1 2*z*w w^2];
Qlpf = tf(num,den);
QFilter = QFilter*ss(Qlpf);

for i = 2:length(filters_q) 

    wn = filters_q(i,1); dn = filters_q(i,2); wd = filters_q(i,3); dd = filters_q(i,4);
    num = [1/wn.^2 2*dn/wn 1];
    den = [1/wd.^2 2*dd/wd 1];
    n = ss(tf(num,den));
    QFilter = QFilter*n;

end

filters_naf = xlsread('Filters.xlsx', 'NAF');

w = filters_naf(1,3);  z= filters_naf(1,4);
num = w^2;
den = [1 2*z*w w^2];
Nlpf = tf(num,den);
NFilter = NFilter*ss(Nlpf);

for i = 2:length(filters_naf) 

    wn = filters_naf(i,1); dn = filters_naf(i,2); wd = filters_naf(i,3); dd = filters_naf(i,4);
    num = [1/wn.^2 2*dn/wn 1];
    den = [1/wd.^2 2*dd/wd 1];
    n = ss(tf(num,den));
    NFilter = NFilter*n;

end

NAF_Gain = 0.045;
NAF = tf(1,[1 0]);
NAF = NAF_Gain*ss(NAF);
NAF = NAF*NFilter;

NAF.inputname = 'Nacelle fore-aft acceleration';
NAF.outputname = 'Collective pitch angle demand';

%% curve fit %%
Order = 1;  %Station Number
S_sign = 1;
FinePitch = PitchAngles(1); %Minimum Fine PitchAngles
PitchLoop = 0; %PitchAngles Loop Control Flag

while ~PitchLoop&&S_sign < length(PitchAngles)
    if PitchAngles(S_sign) > FinePitch
        Q_station = S_sign - 1; %% Torque loop
        Q_pitchangle = PitchAngles(S_sign - 1);
        Q_windspeed = Windspeeds(S_sign - 1);
        
        P_station(Order) = S_sign; %%Pitch loop first station
        P_pitchangle(Order) = PitchAngles(S_sign);
        P_windspeed(Order) = Windspeeds(S_sign);
        
        PitchLoop = 1;
    end
    S_sign = S_sign + 1;
end

Order = Order + 1;
P_station(Order) = S_sign + 1;
P_pitchangle(Order) = PitchAngles(S_sign + 1);  %%Second Station
P_windspeed(Order) = Windspeeds(S_sign + 1);


for j = S_sign + 3:4:length(PitchAngles)
    Order = Order + 1;
    if j < (length(PitchAngles)-4)
        P_station(Order) = j;
        P_pitchangle(Order) = PitchAngles(j);  %%Further Station
        P_windspeed(Order) = Windspeeds(j);
    else
        P_station(Order) = length(PitchAngles);
        P_pitchangle(Order) = PitchAngles(end);  %%Further Station
        P_windspeed(Order) = Windspeeds(end);
    end
end


c_s = 0;
for i = 1:length(PitchAngles)
    if PitchAngles(i)>PitchAngles(1)
        c_s = i;
        break;
    end
end
pitch_angles = PitchAngles(c_s:end);
pitch_angles_dot = zeros(1,length(pitch_angles));
for i = 1:length(pitch_angles)
    if i == 1
        pitch_angles_dot(i) = pitch_angles(i)-PitchAngles(1);
    else
        pitch_angles_dot(i) = pitch_angles(i)-pitch_angles(i-1);
    end
end

P_pitchangle_dot = zeros(1,length(P_pitchangle));
for i  = 1:length(P_pitchangle)
    ind = find(pitch_angles == P_pitchangle(i));
    P_pitchangle_dot(i) = pitch_angles_dot(ind);
end

kp_list = zeros(1,length(P_pitchangle));
ti_list = zeros(1,length(P_pitchangle));
kd_list = zeros(1,length(P_pitchangle));
ts_list = zeros(1,length(P_pitchangle));
tr_list = zeros(1,length(P_pitchangle));

kp_list(1) = first_kp;
ti_list(1) = first_ti;
kd_list(1) = first_kd;
ts_list(1) = first_ts;
tr_list(1) = first_tr;

cal_set = textread('config.txt');
fit_coe = cal_set(1,3);

for i  = 2:length(P_pitchangle)
    coe = power(P_pitchangle_dot(i)/P_pitchangle_dot(1),fit_coe);
    kp_list(i) = coe*kp_list(1);
    ti_list(i) = coe*ti_list(1);
    kd_list(i) = coe*kd_list(1);
    
    Kp = kp_list(i);
    Ti = ti_list(i);
    Kd = kd_list(i);
    Td = first_td;
    
    PPID = ss(tf([Kd + Kp*Td, Kp + Kp/Ti*Td, Kp/Ti], [Td, 1, 0]));
    PPID = PPID*PFilter; 
    PPID.inputname = 'Measured generator speed';
    PPID.outputname ='Collective pitch angle demand';

    % Assign Loop %
    Station = P_station(i);
    iAzimuth = 1;
    sys = ss(SYSTURB.A(:,:,Station,iAzimuth),SYSTURB.B(:,:,Station,iAzimuth),SYSTURB.C(:,:,Station,iAzimuth),SYSTURB.D(:,:,Station,iAzimuth));
    sys.inputname = SYSTURB.inputname;  %A list of all the inputnames
    sys.outputname = SYSTURB.outputname; %A 

    sys = addController(sys, NAF, NAF_type);
    Controller = PPID;
    ControllerType = PPID_type;

    sysd = convertSsTime(sys,ControllerType.TimeStep,'ZOH','Plant');
    sysdol = getSiso(sysd,Controller.outputname,Controller.inputname);
    sysdol.outputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep);
    sysdol.inputdelay = round((ControllerType.InputDelay + ControllerType.OutputDelay)/ControllerType.TimeStep)-sysdol.outputdelay;
    Controller = c2d(Controller,ControllerType.TimeStep,'tustin');   

    Controller.inputname = PPID.inputname;
    Controller.outputname = PPID.outputname;

    ControllerDelay = Controller;
    ControllerDelay.inputdelay = round(ControllerType.InputDelay/ControllerType.TimeStep); % add delays to the controller
    ControllerDelay.outputdelay = round((ControllerType.InputDelay + ControllerType.OutputDelay)/ControllerType.TimeStep)-sysdol.outputdelay; % add delays to the controller

    % The closed loop system %
    syscl = addController(sys, ControllerDelay, ControllerType);
    [y_sam,t_sam] = step(getSiso(syscl,'Collective wind speed','Blade 1 pitch angle'),100);
    Step1 = stepinfo(y_sam,t_sam);
    
    result = [Kp,Ti,Kd,Td,Step1.SettlingTime,Step1.RiseTime,Step1.Overshoot];
    
    ts_list(i) = Step1.SettlingTime;
    tr_list(i) = Step1.RiseTime;

    fid = fopen(['..\Pitch', num2str(i),'\Result.txt'],'wt');
    for k = 1:size(result,1)
        for j = 1:size(result,2)
            fprintf(fid, '%.3f ', result(k,j));
        end
        if i < size(result,1)
            fprintf(fid, '\n');
        end
    end
    fclose(fid);  
end

h = figure('Visible', 'off');
plot(P_pitchangle,ts_list);
xlabel Pitchangle
ylabel SettlingTime
grid on
saveas(h, 'SettlingTime','png')

g = figure('Visible', 'off');
plot(P_pitchangle,tr_list);
xlabel Pitchangle
ylabel RiseTime
grid on
saveas(g, 'RiseTime','png')





    








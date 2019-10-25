clc
close all
clear all
warning off

global PPID_type QPID_type PFilter QFilter

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

%% Load Station %%

% Temp_p = xlsread('..\Temp\Station.xlsx', 'Pitch');
% Station_p = Temp_p(:,1);
% Temp_q = xlsread('..\Temp\Station.xlsx', 'Torque');
% Station_q = Temp_q(1,1);
% Station = Station_p(1);

Info_s = textread('Station.txt');
is_pitch = Info_s(1,1);
Station = Info_s(1,2);
is_first = Info_s(1,3);
p1_td = Info_s(1,4);

load('linmod1.mat');   
iAzimuth = 1;  

sys = ss(SYSTURB.A(:,:,Station,iAzimuth),SYSTURB.B(:,:,Station,iAzimuth),SYSTURB.C(:,:,Station,iAzimuth),SYSTURB.D(:,:,Station,iAzimuth));
sys.inputname = SYSTURB.inputname;  %A list of all the inputnames
sys.outputname = SYSTURB.outputname; %A 

if is_pitch
    sys = addController(sys, NAF, NAF_type);
end

result = GAMultiObj(sys,is_pitch,is_first,p1_td);

if exist('Result.txt','file')
    delete('Result.txt')
end

fid = fopen('Result.txt','wt');
for i = 1:size(result,1)
    for j = 1:size(result,2)
        fprintf(fid, '%.3f ', result(i,j));
    end
    if i < size(result,1)
        fprintf(fid, '\n');
    end
end
fclose(fid);




% April 2013
% Driver for the abrupt CL program.
% Wilhelms did a 30 s pacing and measured alternans of the last 5 oscillations.
% They measured APD50.
% they also had AP profiles from various high and low pacing rates.
% the range of CL they took was 0.2 s to 1 s.
% In out case, we now want to run our program for:
% APD90, APD50, APD30, min cai, max cai.
% pacing for 30s at given CL, not a given number of beats.
% then do measurements for last 10 beats.
% ODE15S works better than ODE23 in this case.
%
% Reduce the program,
% then put in abrupt CL.
% this takes only a few seconds on the local interpreted machine - dont bother with Zen.

clear all;
clear all;
close all;
close all;

% parallel processing does not support global variables
global Istim stim;
global gtomult;
global gk1mult;

global AF;

% global i_to i_K1 i_Ca_L i_Kr i_Ks i_Na;

% set up your simulation here.
drawScrn = 0;

tic

% Lau et al. Cardiovasc Res. 1988, 22: 67-22.
% the baseline CL are (CL1 = CL3): 1200, 545.
% the abrupt changed CL for 1200: 857, 666, 545 (expt1)
% the abrupt changed CL for 545: 1200, 857, 666 (expt2)
% at each pacing rate, the pacing is for 10 minutes.

for CL2_count=3:1:6
for BB=1:2
for AF_it=1:2
for cellType=1:2 % 1 = CRN, 2 = Grandi, 3 = KT


AF = AF_it - 1;

if CL2_count==1
CL1 = 1200.0;
CL2 = 857;
CL3 = CL1;
end;

if CL2_count==2
CL1 = 1200.0;
CL2 = 666;
CL3 = CL1;
end;

if CL2_count==3
CL1 = 1200.0;
CL2 = 545;
CL3 = CL1;
end;

if CL2_count==4
CL1 = 545.;
CL2 = 1200;
CL3 = CL1;
end;

if CL2_count==5
CL1 = 545;
CL2 = 857;
CL3 = CL1;
end;

if CL2_count==6
CL1 = 545.0;
CL2 = 666;
CL3 = CL1;
end;


if BB==1
gtomult = 1.0;
gk1mult = 1.0;
end;

if BB==2
gtomult = 1.0-0.41;
gk1mult = 1.0-0.34;
end;



i_s = 1; % index for the stimulus number
% Ncond in this program is the total number of stimuli applied by CL_i
Ncond = 1000; % tenminutes/CL1;

%% so here we set up the different setting for cell model.
if cellType==1 % CRN
y0 = [0.0, 1.0, 1.0, 1.367e-4, 7.755e-1, 9.996e-1, 9.649e-1, 9.775e-1, ...
    2.908e-3, 1.013e-4, 1.488, 1.488, 1.39e2, 1.117e1, -81.18, 0.0, 1.869e-2, ...
    3.043e-2, 9.992e-1, 4.966e-3, 9.986e-1, 0, 0];
    stim_current = -2000;                % stimulation current
    tstim = 3;               % stimulation duration (tstim > dt)
    oi_V = 15;
    oi_Cai = 10;
    tUnits = 1;
end;

if cellType==2 % Grandi
p = 0;
% dy_grandi improved for no-INaL, no-mc
    y0 = [1.4056270e-03   9.8670050e-01   9.9156200e-01   7.1756620e-06   1.0006810e+00   2.4219910e-02 ...
    1.4526050e-02   4.0515740e-03   9.9455110e-01   4.0515740e-03   9.9455110e-01   8.6413860e-03   ...
    5.4120340e-03   8.8843320e-01   8.1566280e-07   1.0242740e-07   3.5398920e+00   7.7208540e-01   ...
    8.7731910e-03   1.0782830e-01   1.5240020e-02   2.9119160e-04   1.2987540e-03   1.3819820e-01   ...
    2.1431650e-03   9.5663550e-03   1.1103630e-01   7.3478880e-03   7.2973780e-02  1.2429880e+00   ...
    1.0000000e-02   9.1360000e+00   9.1360000e+00   9.1360000e+00   1.2000000e+02   1.7374750e-04   ...
    1.0318120e-04   8.5974010e-05  -7.500000e+01   9.9460000e-01   1.0000000e+00 ...
    0.0000000e+00 1.0000000e+00 0 0];
% size(y0)
    stim_current = 12.5;                % stimulation current
    tstim = 5;               % stimulation duration (tstim > dt)
    oi_V = 39;
    oi_Cai = 37; % this may not be Cai, but in ERP we do not care as of now.
    tUnits = 1;
end;

if cellType==3 % KT
    y0 = [-7.702786e+001 2.775812e-003 9.039100e-001 9.039673e-001 1.060917e-005 9.988566e-001 ...
        9.988624e-001 9.744374e-001 9.594258e-004 9.543380e-001 3.111703e-004 9.751094e-001  ...
        4.109751e-003 4.189417e-005 5.620665e-002 3.975097e-005 9.999717e-001 2.455297e-001 ...
        9.478514e-005 9.993722e-001 1.925362e-001 7.765503e-005 9.995086e-001 2.010345e-001 ... 
        5.674947e-005 9.995604e-001 2.163122e-001 4.638565e-003 4.512078e-003 4.326409e-003 ... 
        4.250445e-003 8.691504e+000 9.286860e+000 1.346313e+002 1.619377e-004 1.354965e-004 ... 
        1.381421e-004 1.442087e-004 1.561844e-004 6.189225e-001 6.076289e-001 5.905266e-001 5.738108e-001];
    stim_current = -80;                % stimulation current
    tstim = 0.030;               % stimulation duration (must be > dt)
    oi_V = 1;
    oi_Cai = 35;
    tUnits = 1e-3;
end;

Istim = stim_current;                % stimulation current

CL1 = CL1*tUnits;
CL2 = CL2*tUnits;
CL3 = CL3*tUnits;
%
% the metadata are set up, now do the conditioning pulses.
%
stim = zeros(Ncond,2);
for i = 0:Ncond
   stim(2*i+1:2*i+2,:) = [CL1*i CL1*i+tstim; Istim 0]';
   timestim(i_s) = CL1*i;
   i_s = i_s + 1;
end;

% obtain solution of conditioning pulses
tspan = [0 CL1*Ncond];
dt = 1.0*tUnits;    % as we are doing conditioning where we do not care for the output as such, we take a large time step. the ode solver is adaptive time stepping.
options = odeset('RelTol',1e-3,'MaxStep',dt);

if cellType==1
[t,y] = ode15s(@dy_crn,tspan,y0,options,'cond');
end;

if cellType==2
[t,y] = ode15s(@dy_grandi,tspan,y0,options,'cond');
end;

if cellType==3
[t,y] = ode15s(@dy_kt,tspan,y0,options,'cond');
end;

for i = 0:Ncond
   stim(2*i+1:2*i+2,:) = [CL2*i CL2*i+tstim; Istim 0]';
if i>0
   timestim(i_s) = CL2*i+CL1*Ncond;
   i_s = i_s + 1;
end;
end;

y1=y;
t1=t;
y0 = y(length(t),:);
% Ncond = 10; % tenminutes/CL2;
% obtain solution of conditioning pulses
tspan = [0 CL2*Ncond];
dt = 1.0*tUnits;    % as we are doing conditioning where we do not care for the output as such, we take a large time step. the ode solver is adaptive time stepping.
options = odeset('RelTol',1e-3,'MaxStep',dt);

if cellType==1
[t,y] = ode15s(@dy_crn,tspan,y0,options,'cond');
end;

if cellType==2
[t,y] = ode15s(@dy_grandi,tspan,y0,options,'cond');
end;

if cellType==3
[t,y] = ode15s(@dy_kt,tspan,y0,options,'cond');
end;
t=t+CL1*Ncond;
y2=y;
t2=t;
y0 = y(length(t),:);

for i = 0:Ncond
   stim(2*i+1:2*i+2,:) = [CL3*i CL3*i+tstim; Istim 0]';
if i>0
   timestim(i_s) = CL3*i+CL1*Ncond+CL2*Ncond;
   i_s = i_s + 1;
end;
end;


% obtain solution of conditioning pulses
% Ncond = tenminutes/CL3;
tspan = [0 CL3*Ncond];
dt = 1.0*tUnits;    % as we are doing conditioning where we do not care for the output as such, we take a large time step. the ode solver is adaptive time stepping.
options = odeset('RelTol',1e-3,'MaxStep',dt);

if cellType==1
[t,y] = ode15s(@dy_crn,tspan,y0,options,'cond');
end;

if cellType==2
[t,y] = ode15s(@dy_grandi,tspan,y0,options,'cond');
end;

if cellType==3
[t,y] = ode15s(@dy_kt,tspan,y0,options,'cond');
end;
t=t+CL1*Ncond+CL2*Ncond;
y3=y;
t3=t;

clear y;
clear t;
y = [y1;y2;y3];
t = [t1;t2;t3];

for i = 1:length(timestim)
    istim(i) = min(find(t>=timestim(i)));
end;

% now measurements,
% now extract the peaks first.
if cellType==1
[overshoot oslocs]=findpeaks( y(:,oi_V),'minpeakheight',10);
[resting restlocs]=findpeaks(-y(:,oi_V),'minpeakheight',30);
resting=-resting;
end;

if cellType==2
[overshoot oslocs]=findpeaks( y(:,oi_V),'minpeakheight',-10);
[resting restlocs]=findpeaks(-y(:,oi_V),'minpeakheight',0);
resting=-resting;
end;

if cellType==3
[overshoot oslocs]=findpeaks( y(:,oi_V),'minpeakheight',-10);
[resting restlocs]=findpeaks(-y(:,oi_V),'minpeakheight',0);
resting=-resting;
end;

[maxCai manCailocs]=findpeaks( y(:,oi_Cai));
[minCai minCailocs]=findpeaks(-y(:,oi_Cai));
minCai=-minCai;

v30=resting*0.3;
v50=resting*0.5;
v90=resting*0.9;

t30 = zeros(1,length(v30)) - 10000.0;
t50 = zeros(1,length(v50)) - 10000.0;
t90 = zeros(1,length(v90)) - 10000.0;

size(t30)
size(v30)
size(resting)

for resc=1:1:length(restlocs)-1 % start at t = restlocs(1).

	% now choose 1 interval between 2 resting times.t(restlocs(resc)) to t(restlocs(resc+1))
	for time_i=restlocs(resc):1:restlocs(resc+1)
		if t30(resc)<0&&y(time_i-1,oi_V)>=v30(resc)&&y(time_i,oi_V)<v30(resc)&&resc<=length(v30)
			t30(resc) = t(time_i); % do this once between every 2 resting times. the last one does not record.
		end;

		if t50(resc)<0&&y(time_i-1,oi_V)>=v50(resc)&&y(time_i,oi_V)<v50(resc)&&resc<=length(v50)
			t50(resc) = t(time_i); % do this once between every 2 resting times. the last one does not record.
		end;

		if t90(resc)<0&&y(time_i-1,oi_V)>=v90(resc)&&y(time_i,oi_V)<v90(resc)&&resc<=length(v90)
			t90(resc) = t(time_i); % do this once between every 2 resting times. the last one does not record.
		end;
	end;
end;

if drawScrn==1
figure;
plot(t,y(:,oi_V)); 
hold on;
plot(t(restlocs),resting,'*b');
hold on;
plot(t(oslocs),overshoot,'+g');
hold on;
plot(t30,v30,'o');

plot(t(istim),y(istim,oi_V),'*r')

axis([0 t(end) -90 50.0]);
drawnow;
end;

apd30=t30(1:end-2)-t(istim(2:length(t30)-1))';
apd50=t50(1:end-2)-t(istim(2:length(t30)-1))';
apd90=t90(1:end-2)-t(istim(2:length(t30)-1))';

apd30_file=apd30(1:end-1)/tUnits;
apd50_file=apd50(1:end-1)/tUnits;
apd90_file=apd90(1:end-1)/tUnits;
maxCai_file=maxCai(1:length(apd30_file));
minCai_file=minCai(1:length(apd30_file));

% apchar=[apd30_file' apd50_file' apd90_file' maxCai_file' minCai_file'];
apchar=[apd30_file' apd50_file' apd90_file'];
restfile =sprintf('abrupt_%d_%g_%g_%d_%d.dat',cellType,CL1,CL2,AF,BB);
save(restfile,'apchar','-ascii','-append');

end; % end CL2_count
end; % celltype
end; % end of AF
end; % end of BB

toc


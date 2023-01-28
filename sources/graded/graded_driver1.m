% July 2013. The idea is to do an AF case for CRN and Grandi at least.
% 
% Driver for the alternans program.
% Wilhelms did a 30 s pacing and measured alternans of the last 5 oscillations.
% They measured APD50.
% they also had AP profiles from various high and low pacing rates.
% the range of CL they took was 0.2 s to 1 s.
% In out case, we now want to run our program for:
% APD90, APD50, APD30, min cai, max cai.
% pacing for 30s at given CL, not a given number of beats.
% then do measurements for last 10 beats.
% ODE15S works better than ODE23 in this case.
clear all;
close all;
clear all;
close all;

% parallel processing does not support global variables
global Istim stim;
global gtomult;
global gk1mult;
global AF;

mult_range = [0.1:0.1:3.0];

% set up your simulation here.
drawScrn = 0;
CL = 1000.0;

tic

for cellType=1:1:1 % 1 = CRN, 2 = Grandi, 3 = KT
for AF_it=1:1:2
for mult_iterategto=1:length(mult_range) % make a loop of this, and cellType if you want, and thats that.
for mult_iterategk1=1:length(mult_range) % make a loop of this, and cellType if you want, and thats that.

AF = AF_it-1;

% 2 parameter map of stuff.
gk1mult = mult_range(mult_iterategk1);
gtomult = mult_range(mult_iterategto);

% Conditioning pulses, and 5 final pulses.
Ncond = 30;

%% so here we set up the different setting for cell model.
if cellType==1 % CRN the crn has INaL
y0 = [0.0, 1.0, 1.0, 1.367e-4, 7.755e-1, 9.996e-1, 9.649e-1, 9.775e-1, ...
    2.908e-3, 1.013e-4, 1.488, 1.488, 1.39e2, 1.117e1, -81.18, 0.0, 1.869e-2, ...
    3.043e-2, 9.992e-1, 4.966e-3, 9.986e-1, 0.0, 0.0];
    stim_current = -2000;                % stimulation current
    tstim = 3;               % stimulation duration (tstim > dt)
    oi_V = 15;
    oi_Cai = 10;
    tUnits = 1;
end;

if cellType==2 % Grandi
p = 0;
% dy_grandi improved for no-Markov chain, but there is an INaL now for the AF case.
    y0 = [1.4056270e-03   9.8670050e-01   9.9156200e-01   7.1756620e-06   1.0006810e+00   2.4219910e-02 ...
    1.4526050e-02   4.0515740e-03   9.9455110e-01   4.0515740e-03   9.9455110e-01   8.6413860e-03   ...
    5.4120340e-03   8.8843320e-01   8.1566280e-07   1.0242740e-07   3.5398920e+00   7.7208540e-01   ...
    8.7731910e-03   1.0782830e-01   1.5240020e-02   2.9119160e-04   1.2987540e-03   1.3819820e-01   ...
    2.1431650e-03   9.5663550e-03   1.1103630e-01   7.3478880e-03   7.2973780e-02  1.2429880e+00   ...
    1.0000000e-02   9.1360000e+00   9.1360000e+00   9.1360000e+00   1.2000000e+02   1.7374750e-04   ...
    1.0318120e-04   8.5974010e-05  -7.500000e+01   9.9460000e-01   1.0000000e+00 ...
    0.0000000e+00 1.0000000e+00 0.0 0.0];
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

CL = CL*tUnits;
%
% the metadata are set up, now do the conditioning pulses.
%
stim = zeros(2*Ncond,2);
for i = 0:Ncond - 1
   stim(2*i+1:2*i+2,:) = [CL*i CL*i + tstim; Istim 0]';
end

% obtain solution of conditioning pulses
tspan = [0 CL*Ncond];
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
			cycle(resc) = CL;
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
plot(t(restlocs),resting,'*');
hold on;
plot(t(oslocs),overshoot,'+');
hold on;
plot(t30,v30,'o');
axis([0 t(end) -80 40.0]);
drawnow;
end;

apd30=t30-t(restlocs)';
apd50=t50-t(restlocs)';
apd90=t90-t(restlocs)';

apd30_file=-1;
apd50_file=-1;
apd90_file=-1;
cycle_file=-1;
maxCai_file=-1;
minCai_file=-1;

if length(apd30)>2
apd30_file=apd30(end-2);
end;
if length(apd50)>2
apd50_file=apd50(end-2);
end;
if length(apd90)>2
apd90_file=apd90(end-2);
end;
if length(cycle)>2
cycle_file=cycle(end-2);
end;
if length(maxCai)>2
maxCai_file=maxCai(end-2);
end;
if length(minCai)>2
minCai_file=minCai(end-2);
end;

apchar=[mult_iterategk1 mult_iterategto gk1mult gtomult cycle_file apd30_file apd50_file apd90_file maxCai_file minCai_file]/tUnits;

% your max Cai and min Cai are not being measured, but this is ERP, so do the important things first.
restfile =sprintf('graded_%d_%d.dat',cellType,AF);
save(restfile,'apchar','-ASCII','-append');

end; % end of mult loop.

end; % end of cellType loop.

end; % end of gk1_or_gto

end; % end of AF loop.

toc


% April 2013
% S1S2 driver. Morph of erp_driver.m

% parallel processing does not support global variables
global Istim stim;
global gtomult;
global gk1mult;
global AF;

close all;
tic
% global i_to i_K1 i_Ca_L i_Kr i_Ks i_Na;

% Conditioning pulses
Ncond = 30;
drawScrn = 0;

% we are doing S1S2 at 2 pacing rates - SR (69 bpm) and BB (58 bpm) pacing rates.
% play in the safe zone so that the program does not fail.
CL_vec = [850 1000 1100];
S2_array = [1000:-25:300];

for BB=1:2

for AF_it=1:2
for cellType=1:1:2 % 1 = CRN, 2 = Grandi, 3 = KT. we do not do KT any more.

AF = AF_it - 1;

if cellType==1||cellType==2
	if BB==1
	gtomult =1.0;
	gk1mult = 1.0;
	end;
	if BB==2
	gtomult = 1.0-0.41;
	gk1mult = 1.0-0.34;
	end;
end;

if cellType==3
gtomult =1.0;
gk1mult = 1.0;
end;

for CL_iterat=1:length(CL_vec) % make a loop of this, and cellType if you want, and thats that.

CL = CL_vec(CL_iterat);

%% so here we set up the different setting for cell model.
if cellType==1 % CRN
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
% dy_grandi improved for no-INaL, no-mc
    y0 = [1.4056270e-03   9.8670050e-01   9.9156200e-01   7.1756620e-06   1.0006810e+00   2.4219910e-02 ...
    1.4526050e-02   4.0515740e-03   9.9455110e-01   4.0515740e-03   9.9455110e-01   8.6413860e-03   ...
    5.4120340e-03   8.8843320e-01   8.1566280e-07   1.0242740e-07   3.5398920e+00   7.7208540e-01   ...
    8.7731910e-03   1.0782830e-01   1.5240020e-02   2.9119160e-04   1.2987540e-03   1.3819820e-01   ...
    2.1431650e-03   9.5663550e-03   1.1103630e-01   7.3478880e-03   7.2973780e-02  1.2429880e+00   ...
    1.0000000e-02   9.1360000e+00   9.1360000e+00   9.1360000e+00   1.2000000e+02   1.7374750e-04   ...
    1.0318120e-04   8.5974010e-05  -7.500000e+01   9.9460000e-01   1.0000000e+00 ...
    0.0000000e+00 1.0000000e+00 0.0 0.0];
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
margin=100.0*tUnits;
%
% the metadata are set up, now do the conditioning pulses.
%
stim = zeros(2*Ncond,2);
for i = 0:Ncond - 1
   stim(2*i+1:2*i+2,:) = [CL*i+margin CL*i + margin + tstim; Istim 0]';
end;


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

if drawScrn==1
figure;
plot(t,y(:,oi_V)); 
drawnow;
end;

%
% now the S1-S2 calculation. put the answer into a data file.
% this gets the data from the conditioning pulses and then makes the
% measurements using the 
%
%
T0 = t(end);
y0_conditioned = y(length(t),:);

for S2_count=1:1:length(S2_array)
S2 = S2_array(S2_count)*tUnits;
% start S2 loop here.
y0=y0_conditioned;
for i = 0:Ncond-1
   stim(2*i+1:2*i+2,:) = [margin+CL*i margin+CL*i+tstim; 0 0]';
end;
% S1
stim(1:2,:) = [margin margin+tstim; stim_current 0]';

% S2
stim(3:4,:) = [S2+margin S2+margin+tstim; stim_current 0]';
options = odeset('RelTol',1e-3,'MaxStep',dt);
tspan = [0 S2+1000.0*tUnits];


if cellType==1
[t,y] = ode15s(@dy_crn,tspan,y0,options,'cond');
end;

if cellType==2
[t,y] = ode15s(@dy_grandi,tspan,y0,options,'cond');
end;

if cellType==3
[t,y] = ode15s(@dy_kt,tspan,y0,options,'cond');
end;

if drawScrn==1
figure;
plot(t,y(:,oi_V)); 
drawnow;
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
% [resting restlocs]=findpeaks(-y(:,oi_V),'minpeakheight',0);
% resting=-resting;
		for t_loc=1:1:length(t)-1
			if t(t_loc)<stim(1)&&t(t_loc+1)>stim(1)
				restlocs(1) = t_loc;
			end;

			if t(t_loc)<stim(3)&&t(t_loc+1)>stim(3)
				restlocs(2) = t_loc;
			end;
		end;

	resting(1) = y(restlocs(1),oi_V);
	resting(2) = y(restlocs(2),oi_V);
end;


[maxCai manCailocs]=findpeaks( y(:,oi_Cai));
[minCai minCailocs]=findpeaks(-y(:,oi_Cai));
minCai=-minCai;

v30=resting*0.3;
v50=resting*0.5;
v90=resting*0.9;

t30 = zeros(1,length(v30)) - 1.0;
t50 = zeros(1,length(v50)) - 1.0;
t90 = zeros(1,length(v90)) - 1.0;

apd30 = zeros(1,length(v30)) - 1.0;
apd50 = zeros(1,length(v50)) - 1.0;
apd90 = zeros(1,length(v90)) - 1.0;

if length(restlocs)==2

	% get your first measurement.
	for first_i=restlocs(1)+1:1:restlocs(2)
		if y(first_i,oi_V)<v30(1)&&y(first_i-1,oi_V)>=v30(1)&&t30(1)<0
		t30(1) = t(first_i);
		end;
		if y(first_i,oi_V)<v50(1)&&y(first_i-1,oi_V)>=v50(1)&&t50(1)<0
		t50(1) = t(first_i);
		end;
		if y(first_i,oi_V)<v90(1)&&y(first_i-1,oi_V)>=v90(1)&&t90(1)<0
		t90(1) = t(first_i);
		end;
	end;
	% and the second. the iterator is an integer.
	for first_i=restlocs(2):1:length(t)
		if y(first_i,oi_V)<v30(2)&&y(first_i-1,oi_V)>=v30(2)&&t30(2)<0
		t30(2) = t(first_i);
		end;
		if y(first_i,oi_V)<v50(2)&&y(first_i-1,oi_V)>=v50(2)&&t50(2)<0
		t50(2) = t(first_i);
		end;
		if y(first_i,oi_V)<v90(2)&&y(first_i-1,oi_V)>=v90(2)&&t90(2)<0
		t90(2) = t(first_i);
		end;
	end;

	apd30=t30-t(restlocs)';
	apd50=t50-t(restlocs)';
	apd90=t90-t(restlocs)';

	di30 = t(restlocs(2))-t30(1);
	di50 = t(restlocs(2))-t50(1);
	di90 = t(restlocs(2))-t90(1);

	caiamp=maxCai-minCai;

	% Interprete this data carefully to make sure you draw the good plots.
	apchar = [CL di30 di50 di90 apd30(2) apd50(2) apd90(2) tUnits*minCai(end) tUnits*maxCai(end) tUnits*caiamp(end)]/tUnits;
	restfile =sprintf('s1s2_rest_%d_%g_%d.dat',cellType,CL,AF);
	save(restfile,'apchar','-ascii','-append');

end; % end of if restlocs if.

end; % end of S2 array.
end; % end of CL loop. end S2 loop here.
end; % end of cellType loop.

end; % end of AF

end; % end of BB

toc


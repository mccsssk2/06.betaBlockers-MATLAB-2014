% April 2013
% erp driver for 3 cell types. CRN, Grandi, KT.
% This program is ready for testing on Zen/other machines.
% this program is for compiling and running on a large machine.
%
% The idea this time is to improve our matlab codes in making them efficient.
% to do this, we use MATLAB functions properly rather than just translate from basic C to MATLAB, and also improved methods.
%
% So this should be integrated driver to run all three cell models. The
% time in KT model is in s while in grandi and CRN in ms. Therefore it uses
% variable tUnits for scaling the time to ms. Now it uses bisection method 
% for the calculation of erp.
% this program can be improved, but I am just going to live with this for the moment.
% the control and BB cases are done.

% parallel processing does not support global variables
global Istim stim;
global gtomult;
global gk1mult;

% global i_to i_K1 i_Ca_L i_Kr i_Ks i_Na;

% Conditioning pulses
Ncond = 20; % used to be 100, but need to hasten it up a bit.
lineProfiles = 0;
drawScrn = 0;

for cellType=1:1:3 % 1 = CRN, 2 = Grandi, 3 = KT

if cellType==1
cl1=[1000:-50:500];
cl2=[490:-5:250];
% set up your simulation here.
gtomult =1.0;
gk1mult = 1.0;
end;
if cellType==2
cl1=[1000:-50:800];
cl2=[790:-10:550];
% set up your simulation here.
gtomult =1.0;
gk1mult = 1.0;
end;
if cellType==3
cl1=[1000:-50:400];
cl2=[390:-10:150];
% set up your simulation here.
gtomult =1.0;
gk1mult = 1.0;
end;

CL_vec = [cl1 cl2];

for CL_iterat=1:length(CL_vec) % make a loop of this, and cellType if you want, and thats that.

CL = CL_vec(CL_iterat);

%% so here we set up the different setting for cell model.
if cellType==1 % CRN
y0 = [0.0, 1.0, 1.0, 1.367e-4, 7.755e-1, 9.996e-1, 9.649e-1, 9.775e-1, ...
    2.908e-3, 1.013e-4, 1.488, 1.488, 1.39e2, 1.117e1, -81.18, 0.0, 1.869e-2, ...
    3.043e-2, 9.992e-1, 4.966e-3, 9.986e-1];
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
    0.0000000e+00 1.0000000e+00];
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

% output file names. dont want this, just pass the array.
% condfile =sprintf('%d_y0__CL%d.dat',cellType,CL);


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
options = odeset('RelTol',1e-4,'MaxStep',dt);

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

% prepare and save output data
% % T0 = t(end)/tUnits;
% % Y0 = y(length(t),:);
% % output = [T0 Y0];
% % save(condfile,'output','-ascii');
% % msg = sprintf('%s saved \n',condfile);
% % disp(msg);
%
%
%
% now the ERP calculation. put the answer into a data file.
% this gets the data from the conditioning pulses and then makes the
% measurements using the 
%
%
T0 = t(end);
Y0 = y(length(t),:);

%% s1s2
% Get the solution. I have reduced the MaxStep.
% definition of indexes - this is for the s1s2 measurements
i_V = 1;
i_Cai = i_V + 1;
i_lastend = i_Cai;

Nbeats = 5;

% Preallocation arrays for output
resting   = zeros(Nbeats,1) + 1000.0;
overshoot = zeros(Nbeats,1) - 1000.0;
t_rest    = zeros(Nbeats,1) - 1000.0;

minCai =    zeros(Nbeats,1) + 1000.0;
maxCai =    zeros(Nbeats,1) - 1000.0;

S1 = T0;

a = 0;
b = 2000*tUnits;

count = 1;

dt = 0.1*tUnits;
options = odeset('RelTol',1e-3,'MaxStep',dt/2);

while 1
    
    S2 = S1 + (a+b)/2; % time at which S2 is applied, not the S1-S2 interval.

    tstart = T0; % what if we use just tstart = 0; then the time do not have to be passed at all.
    t = tstart;
    y0 = Y0;
    j = 1; % counter for the measurements.

    lon = round(((S2-S1)+CL)/dt);

    t_sol = zeros(1,lon);
    y_sol = zeros(lon,i_lastend); % for the voltage variable. Need Cai, Nai, Ki also.

overshoot(1) = -10000.0;
overshoot(2) = -10000.0;

    for i = 1:1:lon

        % the tspan will not work for just t + dt. t + 2*dt is the minimum for ode15s.
        tspan = [tstart:dt/2:tstart+dt];

	if cellType==1
	[t,y] = ode15s(@dy_crn,tspan,y0,options);
	end;

	if cellType==2
	[t,y] = ode15s(@dy_grandi,tspan,y0,options);
	end;

	if cellType==3
	[t,y] = ode15s(@dy_kt,tspan,y0,options);
	end;

        % define Istim here. Again, this is C legacy, but if you want you can vectorise Istim only in this manner.
        % you need a scalar Istim so that you can put it into the Thomas Algorithm/Progonka solver.
        %
        if (t(1) > S1 && t(1) < S1+tstim) || (t(1) > S2 && t(1) < S2 + tstim) %stim_time>0 & stim_time < tstim 
            Istim = stim_current;
        else
            Istim = 0.0;
        end;


        y0 = y(3,:);

        % for now, take the voltage, calcium_i and time
        t_sol(i) = t(3);
        y_sol(i,i_V) = y(3,oi_V);
        y_sol(i,i_Cai) = y(3,oi_Cai);

        % -- do your measurement here.
        if i>1

            if overshoot(j) < y_sol(i,i_V) 
               overshoot(j) = y_sol(i,i_V);
            end;

		if (t_sol(i) >= S2 && j==1)
			j = j + 1;
%			t_sol(i)
%			S1
%			S2
		end;

        end;
        % -- end of measurement if

        % you also want to figure out a way of outputting the currents.
        % increment your time by 1 time steps.
        tstart = tstart + dt;
    end % end of time loop.

if lineProfiles==1
    overshoot
end;

    erp = abs(overshoot(2)/overshoot(1));

    if     (erp < 0.795) % negative
        a = (a+b)/2; b = b;
    elseif (erp > 0.805) % positive
        a = a; b = (a+b)/2;
    else %if (erp >= 0.795 && erp <= 0.805)
        msg = sprintf('erp condition satisfied within %d iterations.',count);
	scrnout = [CL erp S2-S1 S2 S1 a b count];
	scrnout
        disp('');
        break;
    end;

	scrnout = [CL erp S2-S1 S2 S1 a b count];
	scrnout
% for Grandi and KT, I do not get any ERP pacing here.
if drawScrn==1
	figure; plot(t_sol,y_sol(:,i_V)); drawnow;
	close(figure);
end;

% stopping criterion.
if count>100
	break;
end;
    
count = count+1;

if lineProfiles==1
	% choose how much to write.
	k = 1;output = [];
	for i = 1:1:lon
		output(k,:) = [t_sol(i) y_sol(i,:)];
		k = k+1;
	end;
	% warning('please check how is the saving of other files done and update uniformly')
	filename = sprintf('ap_%d_%d_%d.dat',cellType,CL,count);
	save(filename,'output','-ASCII');
	end; % end of if lineProfiles.

end; % end of g loop, the erp calculation loop.

% your max Cai and min Cai are not being measured, but this is ERP, so do the important things first.
restfile =sprintf('erp_rest_%d.dat',cellType);
% The S2-S1 is actually the erp. The variable called erp should always be 0.8.
apcharS2 = [CL (S2-S1)/tUnits erp overshoot(1) overshoot(2)];
save(restfile,'apcharS2','-ascii','-append');

% msg = sprintf('%s saved \n',restfile);
% disp(msg);

a = 0; b = 2000*tUnits;

end; % end of CL loop.

end; % end of cellType loop.


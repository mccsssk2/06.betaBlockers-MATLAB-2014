% 1D for CVr

function oned

matlabpool open

% CL loop
% cellType loop

cellType = 1; % :-1:1 % CellType = 1 is CRN, 2 is Grandi, and 3 is KT
% Set up simulation here, and copy these lines to the dy_* functions.
% I could do it with a function handle, but I dont want to work it out right now.
% put these parameters here as well as in the dy_* functions. The values here are to mark the output.
AF = 0; % AF takes values 0 and 1
BB = 1; % BB takes values 1 and 2

if BB==1
gtomult = 1; % for completeness.
gk1mult = 1; % for completeness.
end;
if BB==2
gtomult = 1-0.41; % for completeness.
gk1mult = 1-0.34; % for completeness.
end;

% for control cases, I need CL_basal to be a bit longer than APD90 and a bit shorter.
for CL_basal=600:-5:100

% global total_current;

t_start = cputime;

if_output = 1; % 1 for wave propagation output, 0 if only wanting CV.

% all your modelling options.
stim_current = 0.0; % give the stimulus a totally different name.
num_beats = 2; % 10 may be better, but then 20 is even better - 5 shows ok adaptation.
beat_counter = 1;

num_cells = 15;
p1 = 5;
p2 = num_cells - 5;

if cellType==1
num_odes=23;   % this depends on the cell model. CRN = 24, Grandi = xx, KT = xy.
v_index = 15; % put the voltage index of the cell model here. CRN = 1, Grandi = 39, KT = ??
dt = 0.1;
dx = 0.1; % in mm.
cm = 100.0; % pF
diffusion = 4.95; % what units what? Same value as used by Pandit in 2005.
gam = diffusion*dt/(cm*dx*dx); % some constant that occurs too many times in the simulation.
CL = CL_basal;
cell_model = @dy_crn;
end;

if cellType==2
num_odes=45; 
v_index = 39; % put the voltage index of the cell model here. CRN = 1, Grandi = 39, KT = ??
dt = 0.1; % Grandi: 0.005 ms
dx = 0.1; % in mm.
diffusion = 0.012;
cm = 1.0; % pF/mm2
gam = diffusion*dt/(cm*dx*dx); % some constant that occurs too many times in the simulation.
CL = CL_basal;
cell_model = @dy_grandi;
end;

if cellType==3
num_odes=43;
v_index = 1; % put the voltage index of the cell model here. CRN = 1, Grandi = 39, KT = ??
dt = 0.0001; % the KT is in seconds.
dx = 0.1; % in micro-m. this may need to be in micro-m since he uses micro-m in his cell code.
diffusion = 4.2; % what units what?
cm = 0.05; % nF. all this works out so that diffusion can be the same in all 3 cases
gam = diffusion*dt/(cm*dx*dx); % some constant that occurs too many times in the simulation, both explicit and implicit.
CL = CL_basal/1000.0;
cell_model = @dy_kt;
end;

y    = zeros(num_cells,num_odes); % so remember the way this is done: the first index is for location, the second is for the ODE.
ynew = zeros(num_cells,num_odes); % so remember the way this is done: the first index is for location, the second is for the ODE.

a = zeros(num_cells,1); % work out if these are rows, or columns.
b = zeros(num_cells,1);
c = zeros(num_cells,1);
rhs = zeros(num_cells,1);
% y0 = zeros(num_odes,1); % becayse we are setting y0 = y(i,:) we dont need
% this creation - it is being created at the good time.
solution = zeros(num_cells,1);

% initialise the complete state. Can we get "steady state" initial conditions for this
for i=1:1:num_cells

if cellType==1 % CRN
y(i,:) = [2.35e-112, 1.0, 0.9992, 1.367e-4, 7.755e-1, 9.996e-1, 9.649e-1, 9.775e-1, ...
    2.908e-3, 1.013e-4, 1.488, 1.488, 1.39e2, 1.117e1, -81.1, 3.296e-5, 1.869e-2, ...
    3.043e-2, 9.992e-1, 4.966e-3, 9.986e-1, 0, 0];
end;

if cellType==2 % Grandi
y(i,:) = [1.4056270e-03   9.8670050e-01   9.9156200e-01   7.1756620e-06   1.0006810e+00   2.4219910e-02 ...
    1.4526050e-02   4.0515740e-03   9.9455110e-01   4.0515740e-03   9.9455110e-01   8.6413860e-03   ...
    5.4120340e-03   8.8843320e-01   8.1566280e-07   1.0242740e-07   3.5398920e+00   7.7208540e-01   ...
    8.7731910e-03   1.0782830e-01   1.5240020e-02   2.9119160e-04   1.2987540e-03   1.3819820e-01   ...
    2.1431650e-03   9.5663550e-03   1.1103630e-01   7.3478880e-03   7.2973780e-02  1.2429880e+00   ...
    1.0000000e-02   9.1360000e+00   9.1360000e+00   9.1360000e+00   1.2000000e+02   1.7374750e-04   ...
    1.0318120e-04   8.5974010e-05  -7.500000e+01   9.9460000e-01   1.0000000e+00 ...
    0.0000000e+00 1.0000000e+00 0 0];
end;

if cellType==3 % KT
y(i,:) = [-7.702786e+001 2.775812e-003 9.039100e-001 9.039673e-001 1.060917e-005 9.988566e-001 ...
    9.988624e-001 9.744374e-001 9.594258e-004 9.543380e-001 3.111703e-004 9.751094e-001  ...
    4.109751e-003 4.189417e-005 5.620665e-002 3.975097e-005 9.999717e-001 2.455297e-001 ...
    9.478514e-005 9.993722e-001 1.925362e-001 7.765503e-005 9.995086e-001 2.010345e-001 ... 
    5.674947e-005 9.995604e-001 2.163122e-001 4.638565e-003 4.512078e-003 4.326409e-003 ... 
    4.250445e-003 8.691504e+000 9.286860e+000 1.346313e+002 1.619377e-004 1.354965e-004 ... 
    1.381421e-004 1.442087e-004 1.561844e-004 6.189225e-001 6.076289e-001 5.905266e-001 5.738108e-001];
end;

end;

% set up the a (supra-diagonal), b (diagonal), c (superdiagonal) elements
% of the system matrix.
for i=1:1:num_cells
    a(i) = -gam;
    b(i) = 1.0 + 2.0*gam;
    c(i) = -gam;
end;

% constant diffusion boundary conditions.
a(1) = 0.0;
c(1) = -2.0*gam;
a(num_cells) = -2.0*gam;
c(num_cells) = 0.0;

stim_time = 0.0;
writeCounter = 0;
writeTime = 0.0;

t1 = zeros(1,num_beats) - 1000.0;
t2 = zeros(1,num_beats) - 1000.0;
cv = zeros(1,num_beats) - 1000.0;

% revise the iterator.

total_time = num_beats*CL;

lon = round(total_time/dt)+100; % the longest possible time. Time loop is stopped with a break statement.

t = 0.0;

options = odeset('RelTol',1e-3,'MaxStep',dt);
ycell = zeros(3,num_odes);

for time_i=1:lon

% put in dt choice here, then calculate the gam again

% solve all odes and set up the RHS. All the i calculations are independent of i+1 and i-1.
% so you can parallelise this loop.
parfor i=1:num_cells
% for i=1:1:num_cells
	y0 = y(i,:);
	tspan = t:dt/2:t+dt;
%	[tloc,ycell] = ode15s(cell_model,tspan,y0,options);
	[tloc,ycell] = ode113(cell_model,tspan,y0,options);
	ynew(i,:) = ycell(end,:);
end; % end of ODE solution.

for i=1:1:num_cells

if cellType==1
	if stim_time<7.0 && i<5
		stim_current = 1200.0;
	else
		stim_current = 0.0;
	end;
	rhs(i) = dt*stim_current/cm + ynew(i,v_index);
end;
if cellType==2
	if stim_time<2.0 && i<8
		stim_current = 1.25*12.5;
	else
		stim_current = 0.0;
	end;
	rhs(i) = dt*stim_current + ynew(i,v_index);
end;
if cellType==3
	if stim_time<0.050 && i<8
		stim_current = 3.0*80.0;
	else
		stim_current = 0.0;
	end;
	rhs(i)=ynew(i,v_index) + dt*stim_current/cm ;
end;

end; % end of for.

% now call the TDMA solver and the solution for dt is done.
solution = TDMAsolver(a,b,c,rhs);
% stuff this solution into v_index of each y.
ynew(:,v_index) = solution;

if cellType == 1 || cellType == 2
	writeTime2 = writeTime;
else
	writeTime2 = writeTime*1000;
end;

if writeTime2>=2.0 && if_output==1
	filename = sprintf('oned_%d_%g_%d_%d.dat',cellType,CL,AF,BB);
	dlmwrite(filename,solution,'delimiter',' ','precision','%10.10f','-append');
	writeTime = 0.0;
end;

if y(p1,v_index)<-40.0 && ynew(p1,v_index)>=-40.0 && t1(beat_counter)<0.0 && t2(beat_counter) < 0.0
	t1(beat_counter) = t;
end;

if y(p2,v_index)<-40.0 && ynew(p2,v_index)>=-40.0 && t2(beat_counter)<0.0 && t1(beat_counter)>0.0
	t2(beat_counter) = t;
end;

if t2(beat_counter)>0.0
	cv(beat_counter) = (p2-p1)*dx/(t2(beat_counter)-t1(beat_counter));
		if beat_counter == num_beats
			break;
		end;
	beat_counter = beat_counter + 1;
end;

y  = ynew; % copy new to old.

stim_time = stim_time + dt;
t = t + dt;
writeTime = writeTime + dt;

	if stim_time > CL
		stim_time = 0.0;
	end;

	if t>=total_time
		break;
	end;

end; % end of time loop.

cv_data = [CL cv];
clear filename;
filename = sprintf('cv_%d_%d_%d_%d.dat',cellType,CL,AF,BB);
dlmwrite(filename,cv_data,'delimiter',' ','precision','%10.10f','-append');

t_finish = cputime - t_start;
runtime_data = [CL num_beats t_finish];
clear filename;
filename = sprintf('cpu_%d_%d_%d_%d.dat',cellType,CL,AF,BB);
dlmwrite(filename,runtime_data,'delimiter',' ','precision','%10.10f','-append');

end;

matlabpool close; % end of parallel.

end

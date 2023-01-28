% driver for the 2D solver.
% starting from oned.m in VW/control
% SK 12 July 2011. Revised 1 April 2012 w.r.t. ADI solver.
% worst case: 20 mins for 1 ms (Grandi) so 2 months for 10 s to finish.
function some_output = twod_init(cellType)

close all;

% all your modelling options.
% cellType = 1; % CellType = 1 is CRN, 2 is Grandi, and 3 is KT
CL_basal = 2000.0; % a total of 2 s is more than plenty.
num_beats = 1;
% num_cells = 375;
num_cells = 188;
outputRate = 1.0; % every 10 time steps, i.e. every ms. That makes 10k outputs. you want all this output because the thing is slow.
dx = 0.2; % mm

tempy=zeros(num_cells,num_cells); % dont do this each time.

if cellType==1
	num_odes=21;                       % this depends on the cell model. CRN = 24, Grandi = xx, KT = xy.
	v_index = 15;                      % put the voltage index of the cell model here. CRN = 1, Grandi = 39, KT = ??
	dt = 0.1;
	cm = 100.0;                        % pF
%	diffusion = 10.0;                  
	diffusion = 4.95;                  % what units? Same value as used by Pandit in 2005.
	gam = diffusion*dt/(2.0*cm*dx*dx); % some constant that occurs too many times in the simulation.
	CL = CL_basal;
	dy_model = @dy_crn;
	init = load('oned_crn_phase.dat');
	file_length = length(init);
	fecg = 'crn_ecg.dat';
	lam = 60;
end;
if cellType==2
	num_odes=44;                       % this depends on the cell model. CRN = 24, Grandi = xx, KT = xy.
	v_index = 39;                      % put the voltage index of the cell model here. CRN = 1, Grandi = 39, KT = ??
	dt = 0.1;
%	diffusion = 4*0.023225;
	diffusion = 0.012;
	cm = 1.0;                          % pF/mm2
	gam = diffusion*dt/(2.0*cm*dx*dx); % some constant that occurs too many times in the simulation.
	CL = CL_basal;
	dy_model = @dy_grandi;
	init = load('oned_grandi_phase.dat');
	file_length = length(init);
	fecg = 'grandi_ecg.dat';
	lam = 60.0;
end;
if cellType==3
	num_odes=43;                       % this depends on the cell model. CRN = 21, Grandi = xx, KT = xy.
	v_index = 1;                       % put the voltage index of the cell model here. CRN = 1, Grandi = 39, KT = ??
	dt = 0.0001;                       % the KT is in seconds.
%	diffusion = 8.5;                   
	diffusion = 4.2;                   % what units?
	cm = 0.05;                         % nF. all this works out so that diffusion can be the same in all 3 cases
	gam = diffusion*dt/(2.0*cm*dx*dx); % some constant that occurs too many times in the simulation, both explicit and implicit.
	CL = CL_basal/1000.0;
	dy_model = @dy_kt;
	init = load('oned_kt_phase.dat');
	file_length = length(init);
	outputRate = 0.1/1000.0;
	fecg = 'kt_ecg.dat';
	lam = 60.0;
end;

y        = zeros(num_cells,num_cells,num_odes); % so remember the way this is done: the first index is for location, the second is for the ODE.
ynew     = zeros(num_cells,num_cells,num_odes); % so remember the way this is done: the first index is for location, the second is for the ODE.
%a        = zeros(num_cells,1);                  % work out if these are rows, or columns.
%b        = zeros(num_cells,1);
c        = zeros(num_cells,1);
%rhsloop  = zeros(num_cells,num_cells);
rhs	 = zeros(num_cells,1);
solution = zeros(num_cells,1);
potential= zeros(num_cells,num_cells);
potentialnew = zeros(num_cells,num_cells);
dy = zeros(num_cells,num_cells);

c(1:1:end) = -gam;

c(1)         = -2.0*gam;
% a(num_cells) = -2.0*gam;
c(num_cells) = 0.0;

% now do the phase.

centre_x=num_cells*dx/2.0+0.05;
centre_y=num_cells*dx/2.0+0.05;

phaseu=zeros(num_cells,num_cells);

for i=1:1:num_cells
for j=1:1:num_cells
	x_loc=i*dx; y_loc=j*dx;
	phaseu(i,j)=atan2((y_loc-centre_y),(x_loc-centre_x))-2*pi*hypot(x_loc-centre_x,y_loc-centre_y)/lam;
end;
end;

figure;
sizzee = size(phaseu);
xx=1:1:sizzee(2);
yy=1:1:sizzee(1);
[X,Y]=meshgrid(xx,yy);
surf(X,Y,phaseu);
view(0,90);
% caxis([0.0 2*pi]);
shading interp;
colormap(jet);
grid off;
colorbar;

phaseu = mod(phaseu,2*pi);
max(max(phaseu))
min(min(phaseu))
phaseu = phaseu/max(max(phaseu));
max(max(phaseu))
min(min(phaseu))

figure;
sizzee = size(phaseu);
xx=1:1:sizzee(2);
yy=1:1:sizzee(1);
[X,Y]=meshgrid(xx,yy);
surf(X,Y,phaseu);
view(0,90);
% caxis([0.0 2*pi]);
shading interp;
colormap(jet);
grid off;
colorbar;

for i=1:1:num_cells
for j=1:1:num_cells
	temp_k=uint64(file_length*phaseu(i,j)); % and now the line numbers are non-negative numbers.
	if temp_k<1
		temp_k = 1;
	end;
	if cellType==1 || cellType==3
	y(i,j,:)=init(temp_k,:); % special case for Grandi.
	else
	y(i,j,1:41) = init(temp_k,1:41);
	y(i,j,42:44) = init(temp_k,58:60);
	end;
end;
end;

tempy = y(:,:,v_index);

figure;
sizzee = size(tempy);
xx=1:1:sizzee(2);
yy=1:1:sizzee(1);
[X,Y]=meshgrid(xx,yy);
surf(X,Y,tempy);
view(0,90);
shading interp;
colormap(jet);
grid off;
colorbar;


% driver for the 2D solver.
% starting from oned.m in VW/control
% SK 12 July 2011. Revised 1 April 2012 w.r.t. ADI solver.
% worst case: 20 mins for 1 ms (Grandi) so 2 months for 10 s to finish.
function some_output = twod(cellType)

close all;
if matlabpool('size') > 0
matlabpool close;
end;
% 
matlabpool; % leave it so, that it goes and finds the max number of procs (8 is max). 3 machines available: Volks: 4 procs, Cefalu:24 procs, Winfree: 8 procs

% all your modelling options.
cellType = 2; % CellType = 1 is CRN, 2 is Grandi, and 3 is KT
CL_basal = 5000.0; % a total of 2 s is more than plenty.
num_beats = 1;
% num_cells = 375;
num_cells = 188;
outputRate = 1.0; % every 10 time steps, i.e. every ms. That makes 10k outputs. you want all this output because the thing is slow.
dx = 0.2; % mm, implicit solver, so step can be bigger.

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
	lam = 3000;
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
	lam = 3000.0;
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
	outputRate = 1.0/1000.0;
	fecg = 'kt_ecg.dat';
	lam = 3000.0;
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
phaseu = mod(phaseu,2*pi);
phaseu = phaseu/max(max(phaseu));

for i=1:1:num_cells
for j=1:1:num_cells
	temp_k=uint64(file_length*phaseu(i,j)); % and now the line numbers are non-negative numbers.
	if temp_k<1
		temp_k = 1;
	end;
%	if cellType==1 || cellType==3
	y(i,j,:)=init(temp_k,:); % special case for Grandi.
%	else
%	y(i,j,1:41) = init(temp_k,1:41);
%	y(i,j,42:44) = init(temp_k,58:60);
%	end;
end;
end;

m_minValue = -85;
m_maxValue = 20;
filename = sprintf('twod_init.tif'); % jpeg is a slightly cheaper format than tif in terms of file size.
tempy=zeros(num_cells,num_cells);
tempy=y(:,:,v_index);
tempy = 255.0 * (tempy - m_minValue) / (m_maxValue - m_minValue);
m8 = uint8(tempy);
% Convert to rgb.
rgbImage = ind2rgb(m8,jet);
shading interp;
% Save the RGB image to disk.
imwrite(rgbImage, filename); 
clf

writeTime = 0.0;
imagecounter = 1;

% the time iterator.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_time = (num_beats)*CL;
lon = round(total_time/dt);
t = 0.0;

tspan = [0:dt/2:dt]; % tspan is simply [0 dt/2 dt] i.e. a constant vector
options = odeset('RelTol',1e-3,'MaxStep',dt);

for time_i=1:lon

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%solve ODEs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

parfor j=1:num_cells
for i=1:num_cells
	y0 = y(i,j,:);
	[~,ycell,s] = ode23(dy_model,tspan,y0,options); % do I really need this ode solver. first choice: RK4.
	ynew(i,j,:) = ycell(end,:);
	dy(i,j) = (ycell(end,v_index) - ycell(1,v_index))/(dt); % a form of total_current, and you dont need to store all the derivatives either.
end;
end;

y = ynew;
if max(max(isnan(y))) > 0 break; end;
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PDE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
potential = y(1:num_cells,1:num_cells,v_index);
for i = 1:2
    for j=1:num_cells
                        top=j+1; bottom = j-1;
        if j==1         top=j+1; bottom = j+1; end;
        if j==num_cells top=j-1; bottom = j-1; end;
         rhs=potential(1:num_cells,j)+gam*(potential(1:num_cells,top)+potential(1:num_cells,bottom)-2.0*potential(1:num_cells,j))+(dt/2.0)*dy(1:num_cells,j);
         potentialnew(:,j)=TDMAsolver(gam,c,rhs);
    end; % end of first dt/2 solution.
    potential = transpose(potentialnew);
    dy = transpose(dy);
end;
toc
y(:,:,v_index) = potential;

if max(max(isnan(y))) > 0 break; end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = t + dt;
writeTime = writeTime + dt;

if writeTime > outputRate
tic
	filename = sprintf('twod_%06d.tif',imagecounter);
	tempy=y(:,:,v_index);
	tempy = 255.0 * (tempy - m_minValue) / (m_maxValue - m_minValue);
	m8 = uint8(tempy);
	% Convert to rgb.
	rgbImage = ind2rgb(m8,jet);
	% Save the RGB image to disk.
	imwrite(rgbImage, filename); 
        imwrite(m8,strcat('gray_',filename));
	clf
	imagecounter = imagecounter + 1;

%	f = figure('Visible','off');
%	sizzee = size(y(:,:,v_index));
%	xx=1:1:sizzee(2);
%	yy=1:1:sizzee(1);
%	[X,Y]=meshgrid(xx,yy);
%	surf(X,Y,y(:,:,v_index));
%	view(0,90);
%	caxis([-85 20]);
%	shading interp;
%	colormap(jet);
%	grid off;
%	colorbar;
%	saveas(f,strcat(filename,'.jpg'),'jpg');

	writeTime = 0.0;


% stopping criterion. This stops the thing too soon.
if max(max(y(:,:,v_index))) < -65.0
	break;
end;
toc
end; % end of output

end; % end of time loop.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matlabpool close;
some_output = 1.0;

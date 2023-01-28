% POTENTIAL = IMG2VM(FILENAME) read gray scale image from the file 
% specified by FILENAME and return potential distribution.

function some_output = tiptraces(cellType)

% any one of these 2 variables are true. not both.
if_tips = 1;
if_traces = 0;

% minimum and maximum potential value represented in the images
m_minValue = -85;
m_maxValue = 20;
num_cells = 188;
counter = 1;
delay = 20;
min_distance = 2.0; % square of some arbitrary but small distance that I like. Right now, I like 2 units.
t=20;
diag_counter = 1;
diag7 = 0.0;
tdiag7 = 0.0;
oldtdiag7 = 0.0;

start=1+delay;
increment=1;

% for tips, the program fails at endfile less than given. so do this manually.
% basal CRN
if cellType==1
endfile=629;
path='results/basal/CRN/';
% -50 for CRN and KT, -30 for Grandi.
iso_value = -50.0; % value for the isoline.
end;

% basal Grandi
if cellType==2
endfile=3443;
path='results/basal/grandi/';
% -50 for CRN and KT, -30 for Grandi.
iso_value = -30.0; % value for the isoline.
end;

% basal KT
if cellType==3
endfile=5000;
path='results/basal/KT/';
% -50 for CRN and KT, -30 for Grandi.
iso_value = -50.0; % value for the isoline.
end;

% bb CRN
if cellType==4
endfile=777;
path='results/BB/CRN/';
% -50 for CRN and KT, -30 for Grandi.
iso_value = -50.0; % value for the isoline.
end;

% basal Grandi
if cellType==5
endfile=3909;
path='results/BB/grandi/';
% -50 for CRN and KT, -30 for Grandi.
iso_value = -30.0; % value for the isoline.
end;

% basal KT
if cellType==6
endfile=5000;
path='results/BB/KT/';
% -50 for CRN and KT, -30 for Grandi.
iso_value = -50.0; % value for the isoline.
end;

iso_value_temp = 255.0 * (iso_value - m_minValue) / (m_maxValue - m_minValue);
iso_value_int = uint8(iso_value_temp);

[x y] = meshgrid([1:1:num_cells],[1:1:num_cells]);
C = zeros(num_cells,num_cells);
Cold = zeros(num_cells,num_cells);
dV = zeros(num_cells,num_cells);
dV_contour = zeros(num_cells,num_cells);

for counter=start:increment:endfile
tic
file = sprintf('gray_twod_%06d.tif',counter);
filename = strcat(path,file);
A = imread(filename);
doubleA = double(A);
potential = doubleA*(m_maxValue - m_minValue)/255 + m_minValue;

file2 = sprintf('gray_twod_%06d.tif',counter-delay);
filename2 = strcat(path,file2);
Aold = imread(filename2);

% plot it, just so you know whats in potential. This may also act as alternate production tifs.
% h1 = figure('visible','off');
% surf(x,y,potential);
% view(2);
% axis tight;
% shading interp;
% axis off;
% colorbar;
% thetif = sprintf('toimg%d.tif',counter);
% thewholetif = strcat(path,thetif);
% saveas(h1,thewholetif,'tif');
% close(h1);

if if_traces==1

% first task, save some sample values. you can fft the output file if required.
diag1 = potential(10,10);
diag2 = potential(40,40);
diag3 = potential(60,60);
diag4 = potential(94,94);
diag5 = potential(188-60,188-60);
diag6 = potential(188-40,188-40);
old_diag7 = diag7;
diag7 = potential(188-10,188-10);
traces_output = [t diag1 diag2 diag3 diag4 diag5 diag6 diag7];
datafile = sprintf('traces_%d.dat',cellType);
wdatafile = strcat(path,datafile);
dlmwrite(datafile,traces_output,'-append','delimiter','\t');
if old_diag7 < iso_value & diag7 >= iso_value
oldtdiag7 = tdiag7;
tdiag7 = t;
diag_counter = diag_counter + 1;
crossingfile = sprintf('crossing_%d.dat',cellType);
cross_output = [t tdiag7 tdiag7-oldtdiag7];
dlmwrite(crossingfile,cross_output,'-append','delimiter','\t');
end;

end;

if if_tips==1

%tip tracing
% potential is a derived matrix of A, or doubleA. so work with A.
% A has values between 0 and 255. The iso_value_int is some value between 0 and 255, which may or may not be part of 
%
[L1, V1, O1] = isocontour(A,iso_value_int);
v1size = size(V1);
v1max = v1size(1);
[L2, V2, O2] = isocontour(Aold,iso_value_int);
v2size = size(V2);
v2max = v2size(1);

% h2 = figure; % ('visible','off');
% plot(V1(:,1),V1(:,2),'.',V2(:,1),V2(:,2),'.');
% axis([1 188 1 188]);
% hold on;

% this is choosing vertices that are close to each other. This is because I do not know how to order the vertices.
counter2 = 1;
tX1 = [];
tY1 = [];
tX2 = [];
tY2 = [];
for v1c=1:1:v1max
	for v2c=1:1:v2max
	distance = (V1(v1c,1) - V2(v2c,1))*(V1(v1c,1) - V2(v2c,1))+(V1(v1c,2) - V2(v2c,2))*(V1(v1c,2) - V2(v2c,2));
		if distance < min_distance
			tX1(counter2) = V1(v1c,1);
			tY1(counter2) = V1(v1c,2);
			tX2(counter2) = V2(v2c,1);
			tY2(counter2) = V2(v2c,2);
			counter2 = counter2 + 1;
		end;
	end;
end;


[X0 Y0] = intersections(tX1,tY1,tX2,tY2);
%% h3 = figure;
% plot(X0,Y0,'+');
% axis([1 188 1 188]);
% hold on;

tip_output = [X0 Y0];
tips = sprintf('tips_%d.dat',cellType);
tipf = strcat(path,tips);
dlmwrite(tips,tip_output,'-append','delimiter','\t');
end;

t
t = t + 1;
cellType
toc
end; % end of loop.

some_output = 1;

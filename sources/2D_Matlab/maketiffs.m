% Program for making the 3 or 4 tiffs per case which show the "good" colormap.
for cellType=1:1:6
% minimum and maximum potential value represented in the images
m_minValue = -85;
m_maxValue = 20;
num_cells = 188;

start=100;
increment=100;

% basal CRN
if cellType==1
endfile=600;
path='results/basal/CRN/';
% -50 for CRN and KT, -30 for Grandi.
iso_value = -50.0; % value for the isoline.
end;

% basal Grandi
if cellType==2
endfile=1000;
path='results/basal/grandi/';
% -50 for CRN and KT, -30 for Grandi.
iso_value = -30.0; % value for the isoline.
end;

% basal KT
if cellType==3
endfile=1000;
path='results/basal/KT/';
% -50 for CRN and KT, -30 for Grandi.
iso_value = -50.0; % value for the isoline.
end;

% bb CRN
if cellType==4
endfile=700;
path='results/BB/CRN/';
% -50 for CRN and KT, -30 for Grandi.
iso_value = -50.0; % value for the isoline.
end;

% basal Grandi
if cellType==5
endfile=1000;
path='results/BB/grandi/';
% -50 for CRN and KT, -30 for Grandi.
iso_value = -30.0; % value for the isoline.
end;

% basal KT
if cellType==6
endfile=1000;
path='results/BB/KT/';
% -50 for CRN and KT, -30 for Grandi.
iso_value = -50.0; % value for the isoline.
end;

[x y] = meshgrid([1:1:num_cells],[1:1:num_cells]);

for counter=start:increment:endfile
	tic
	file = sprintf('gray_twod_%06d.tif',counter);
	filename = strcat(path,file);
	A = imread(filename);
	doubleA = double(A);
	potential = doubleA*(m_maxValue - m_minValue)/255 + m_minValue;

	% plot it, just so you know whats in potential. This may also act as alternate production tifs.
	h1 = figure('visible','off');
	surf(x,y,potential);
	view(2);
	axis tight;
	shading interp;
	axis off;
	colorbar;
	thetif = sprintf('prodfig2d_bbwork_%d_%d.tif',cellType,counter);
%	thewholetif = strcat(path,thetif);
%	saveas(h1,thewholetif,'tif');
	saveas(h1,thetif,'tif');
	close(h1);
	toc;
end;
end;

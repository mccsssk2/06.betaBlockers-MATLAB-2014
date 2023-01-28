close all;
close all;
clear all;
clear all;

cellType=3;

if cellType==1
init = load('oned_crn_phase.dat');
file_length = length(init);
file_length
v_index=15;
end;

if cellType==2
init = load('oned_grandi_phase.dat');
file_length = length(init);
file_length
v_index=39;
end;

if cellType==3
init = load('oned_kt_phase.dat');
file_length = length(init);
file_length
v_index=1;
end;

num_cells=375;
dx=0.1;
centre_x=num_cells*dx/2.0+dx/2;
centre_y=num_cells*dx/2.0+dx/2;
lam=5.0;
phaseu=zeros(num_cells,num_cells);

for i=1:1:num_cells
for j=1:1:num_cells
	x_loc=i*dx;
	y_loc=j*dx;
	phaseu(i,j)=atan2((y_loc-centre_y),(x_loc-centre_x))-2*pi*hypot(x_loc-centre_x,y_loc-centre_y)/lam;
end;
end;

%for i=1:1:num_cells
%for j=1:1:num_cells
%	if phaseu(i,j)>=2*pi phaseu(i,j) = phaseu(i,j) - 2*pi; end;
%	if phaseu(i,j)<=0.0   phaseu(i,j) = phaseu(i,j) + 2*pi; end;
%end;
%end;

phaseu = mod(phaseu,2*pi);

phaseu = phaseu/max(max(phaseu));

for i=1:1:num_cells
for j=1:1:num_cells
	temp_k=uint64(file_length*phaseu(i,j)); % and now the line numbers are non-negative numbers.

	if temp_k<1
		temp_k = 1;
	end;

%	if temp_k>file_length
%		temp_k
%		temp_k = file_length;
%	end;

	y(i,j,:)=init(temp_k,:);
end;
end;

min(min(phaseu))
max(max(phaseu))

figure;
plot(phaseu);

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

figure;
sizzee = size(y(:,:,v_index));
xx=1:1:sizzee(2);
yy=1:1:sizzee(1);
[X,Y]=meshgrid(xx,yy);
surf(X,Y,y(:,:,v_index));
view(0,90);
% caxis([-85 20]);
shading interp;
colormap(jet);
grid off;
colorbar;


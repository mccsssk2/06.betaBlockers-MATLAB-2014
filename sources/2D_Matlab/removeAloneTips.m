% Program to clean up tips_*.dat for plotting with SigmaPlot
function someoutput = removeAloneTips(cellType)

filename = sprintf('tips_%d.dat',cellType);

% load the data
tips = load(filename);
len = length(tips);

% a square of distance from the point to be analyzed
d0 = 1;
d1 = 10;
d2 = 25;
d3 = 100;

for ip = 1:len;
        d0_counter = -1;
        d1_counter = -1;
        d2_counter = -1;
        d3_counter = -1;
    for i = 1:len
        distance = (tips(i,1) - tips(ip,1))^2 + (tips(i,2) - tips(ip,2))^2;
        if distance < d0
            d0_counter = d0_counter + 1;
        end;
        if distance < d1
            d1_counter = d1_counter + 1;
        end;
        if distance < d2
            d2_counter = d2_counter + 1;
        end;
        if distance < d3
            d3_counter = d3_counter + 1;
        end;
    end;
    tip_D(ip,:) = [d0_counter d1_counter d2_counter d3_counter];
end;

% remove points
i = 1; % index of distance to consider
n = 1; % more than how many points have to be around

tips_new = tips(tip_D(:,i) > n,:);

figure;
plot(tips_new(:,1),tips_new(:,2),'.');

filename =strcat('new_',filename);
save(filename,'tips_new','-ascii');

someoutput = 1;


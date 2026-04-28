%%
clear all; clc;
f = fopen('./models_with_state.txt');
models = {};cnt=1;
while ~feof(f)
    l = split(fgetl(f));
    models(cnt,1) = l(1); models(cnt,2) = l(2);
    cnt = cnt+1;
end
fclose(f);
%%
loop = 5; % normalized RMSD
cnt = 0;
target_dir = 'CSol_TMD_situ_all/';
dist = zeros(length(models),length(models));
for i=1:length(models)
    v = dir([target_dir '/*' models{i,1} '_*.txt']);
    for j=1:length(v)
        fp = fopen([target_dir '/' v(j).name]);l=split(fgetl(fp)); fclose(fp);
        dist(i,find(contains(models(:,1), v(j).name(6:end-4)))) = abs(str2double(l{loop}));
        dist(i,i)=0;
        dist(find(contains(models(:,1), v(j).name(6:end-4))),i) = abs(str2double(l{loop}));
        cnt = cnt+1;
    end
end
%% nonmetric
D = dist;
dissimilarities = squareform(D, 'tovector');
[Y,stress,disparities] = mdscale(dissimilarities,2);
distances = pdist(Y);
%%
[dum,ord] = sortrows([disparities(:) dissimilarities(:)]);
plot(dissimilarities,distances,'bo', ...
     dissimilarities(ord),disparities(ord),'r.-');
xlabel('Dissimilarities')
ylabel('Distances/Disparities')
legend({'Distances' 'Disparities'}, 'Location','NorthWest');
%% plotting
nr_my = 19;
s = scatter(Y(1:end-nr_my-1,1), Y(1:end-nr_my-1,2), 40, str2double(models(1:end-nr_my-1,2)), 'filled');
text(Y(1:end-nr_my-1,1)+0.05,Y(1:end-nr_my-1,2),models(1:end-nr_my-1,1), 'Color','blue','FontSize',20);
s.SizeData = 1600;s.MarkerEdgeColor = 'w';s.LineWidth = 1;
%c = colorbar('southoutside', 'Ticks', [1,2,3,4], 'TickLabels', {'Closed','Primed','Open', 'Inactivated'});c.FontSize = 40;
%title("CSol/TMD", 'FontSize',40)

%models{end-9,1} = 'primed_2' 
%models{end-8,1} = 'primed'
%models{end-7,1} = 'open'
%models{end-6,1} = 'apo';
%models{end-5,1} = 'open_ry';
%models{end-4,1} = 'primed_1';
%models{end-3,1} = 'open_2';
%models{end-2,1} = 'apo_1';
%models{end-1,1} = 'apo_2';
%models{end,1}   = 'apo_3';

hold on;
s2 = scatter(Y(end-nr_my:end,1),Y(end-nr_my:end,2), 40, 'r', 'filled');
text(Y(end-nr_my:end,1)+0.05,Y(end-nr_my:end,2), models(end-nr_my:end,1),'Color','blue','FontSize',20);
s2.SizeData = 1600;s2.MarkerEdgeColor = 'w';s2.LineWidth = 1;

%colormap(jet(4));

%axis('equal');
%%
writematrix(Y,'./CSol_TMD.csv');
%%
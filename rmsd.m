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
%%
writematrix(Y,'./CSol_TMD.csv');
%%

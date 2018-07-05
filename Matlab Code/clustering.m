%% clustering
clear
load data filename names
load([filename '_tnorm.mat']);
delete([filename '_cluster_measures.csv'],[filename,'_cluster_vecs.csv'])
if strcmp(variables_to_include,'auto')==1 || isempty(variables_to_include)==1
    variables=find(mean(cdf)<pvalue)+1;
else
    variables=sfind(tnorm(1,:),variables_to_include');
end

t2=cell2mat(tnorm(2:end,variables));
num_of_var=number_of_independent_varibles;
t3=reshape(t2,size(t2,1)/num_of_var,size(t2,2)*num_of_var);
r=1
methods={'neural','kmeans','average','centroid','complete','median','single','weighted','ward'};
tic
r=0;
for j=1:size(methods,2)
    if strcmp(methods{j},'kmeans')==1
        metrics={'sqeuclidean','cityblock','cosine','correlation'};
    elseif strcmp(methods{j},'neural')==1
        metrics={'network som'};
    else
        metrics={'euclidean','seuclidean','cityblock','minkowski','chebychev','cosine','spearman'};
    end
    for jj=1:size(metrics,2)
        %try
        z=1;%greedy_subspace(t2,notes,metrics{jj},methods{j},filename2);
        [~,~,~,~,~,~,~,f,~,fsimilar,T]=dendqual5(t3,z,metrics{jj},methods{j},r,names,[filename,'_cluster_measures.csv']);
        r=1+r
        dlmwrite([filename '_cluster_vecs.csv'],[f,fsimilar,T'],'delimiter', ',' ,'-append')
        clusters(r,:)=[f,fsimilar,T'];
        %end
    end
end
save clusters.mat clusters names t3 f fsimilar filename
toc
'average fscore='
mean(clusters(:,1:2))

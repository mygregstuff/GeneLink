function [v,I,C,d,gd,M,m,f,sgd,fsimilar,T]=dendqual5(t,x,metric,method,r,strains,filename)
[H,PE,norm,ind,nmi,H_clust,H_class,w,c,PC,sgv,v,I,C,d,gd,M,m,f,sgd,fsimilar]=deal(0);

t = x.*t;
t(isnan(t))=0;
notempty = cellfun(@any, strains);
t=t(notempty,:);
genes = strains(notempty)';
u = unique(genes);
num_of_clusters=length(u);


if strcmp(method,'kmeans')==1
    [T,~,~,D]=kmeans(t,length(u),'Distance',metric,'Replicates',10);
elseif strcmp(method,'neural')==1
    net = selforgmap([5 5]);
    net = closeloop(net);
    [net,tr] = train(net,t');
    y = net(t');
    T = vec2ind(y)';
else
    D = pdist(t,metric);
    D(isnan(D))=max(D);
    tree = linkage(D,method);
    outperm = optimalleaforder(tree,D);
    %outperm=1:sum(notempty);
    M=squareform(D);
    M=M(outperm,outperm);
    C = cophenet(tree,D);
    d = abs(diff(M));
    d = mean(d(:));
    %genes = genes(outperm);
    u = unique(genes);
    I = mean(inconsistent(tree,num_of_clusters));
    I = I(1,4);
    T = cluster(tree,num_of_clusters);
end
vec=zeros(1,length(genes));
for k=1:length(u)
    ind=find(ismember(genes,u{k}));
    vec(ind)=k;
end
uvec=1:length(u);
%{
accuracy = (TP + TN) ./ (TN + TP + FN + FP);
sensitivity = TP ./ (TP + FN);
specificity = TN ./ (FP + TN);
precision = TP ./ (TP + FP);
f_score = nansum(nansum( 2*TP ./ (2 .* TP + FP + FN)))/length(vec);
f_score = nansum(norm' .* max( 2*TP ./ (2*TP + FP + FN) ,[],2));
Jaccardind = nansum(norm' .* max( TP ./ (TP + FP + FN) ,[],2));
Fowlkesmallows = nansum( norm' .* max( TP.*sqrt( 1 ./ (TP + FP)./( TP + FN ) ) ,[],2));
nc2 = nchoosek(length(vec),2);
adjrandind =  nansum( norm' .* max( (nc2.*( TP + TN ) - ( TP + FP ).*( TP + FN ) - ( FN + TN ).*( FP + TN )) ./ (nc2^2 - ( TP + FP ).*( TP + FN ) - ( FN + TN ).*( FP + TN )),[],2));
PE = nansum(nansum(PE));
PC = nansum(nansum(PC));
%}

f = fscore(T,vec);
fsimilar=0;
try
    fsimilar = fscore_similar2(genes,T,vec);
end
if (strcmp(method,'kmeans')==0 && strcmp(method,'neural')==0)
    % 'hi'
    % end
    sgd = similargenedist(M,genes');
    [gd,v]=genedist(M(outperm,outperm),vec(outperm));
    %     info=nansum(nansum(TP./length(vec).*log(length(vec).*TP./w./c)));
    %     nmi=2*info/(H_class+H_clust);
    m=[]; m2=[];
    for i=uvec
        ind=find(vec==i);
        m2(:,i)=mean(M(:,ind),2);
        %    outpermm{i}=u{i};
    end
    for i=uvec
        ind=find(vec==i);
        m(i,:)=mean(m2(ind,:),1);
    end
    m3=m;
    m3(uvec,uvec)=0;
end
if exist('filename')>0
    result={metric,method,'f_score',f,'f_score_similar',fsimilar,'Avg_gene_dist',gd,'simmialr gene dist',sgd,'variance',v,'inconsistency',I,'cophenet',C,'roughness',d,num_of_clusters,r};
    fileID = fopen(filename,'a');
    fmt ='%s, %s, %s, %g, %s, %g, %s, %g, %s, %g, %s, %g, %s, %g, %s, %g, %s, %g, %%, %g, %g, %% \r\n';
    fprintf(fileID,fmt, result{:});
    fclose(fileID);
end

end
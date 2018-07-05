clear all;
run Settings.m
load clusters.mat


%create cluster vecs
u = unique(names);
vec=zeros(1,length(names));
for k=1:length(u)
    ind=find(contains(names,u{k}));
    vec(ind)=k;
end
a=double(jet);
ind=round(vec/length(u)*64);
col=a(ind,:);



%%
c=clusters;
A=zeros(size(c,2)-2);
algos=size(c,1); %number of cluster algorithms
for jjj=1:algos
    c3=c(jjj,3:end);
    c1=c(jjj,1);
    c2=c(jjj,2);    
    for j=1:size(c,2)-2
        for jj=1:size(c,2)-2
            %create adjacency matrix
            A(j,jj)=(c1*c2)^sensitivity*double(isequal(c3(j),c3(jj)))+A(j,jj);
        end
        A(j,j)=0;
    end    
end
A=A/max(A(:));
D=A;
%%
figure;
%create color scale
j=linspace(prctile(A(A>.85),0),1,20);
col2=winter;
nodenames=[];
for i=1:length(names)
    ii=num2str(i);
    nodenames{i}=strcat(names{i},'_',ii);
end
k=0;
A2=A
for i=1:length(j)
    A2(A2<j(i))=0;
    G=graph(A2,nodenames);
    if k==0
        h=plot(G,'NodeLabel',G.Nodes.Name,'LineWidth',.5); hold on
        k=1
    end
    highlight(h,G,'EdgeColor',col2(round(i*63/20)+1,:),'LineWidth',4)
    
    %   A(A<.50)=0;
    %   G3=graph(A,names1);
    %   highlight(h,G3,'EdgeColor','r','LineWidth',2)
end

highlight(h,1:size(vec,2),'MarkerSize',7);

for i=1:length(vec)
    highlight(h, i,'NodeColor',col(i,:));
    %   [T,p] = minspantree(G,'Root',i);
    %   highlight(h,T,'EdgeColor','r','LineWidth',1,'LineStyle','--')
end
colormap(winter);
h=colorbar;
ylabel(h,'Relative number of Links Found in Clusters')
set(gcf,'Position',[500 300 700 700]);
%  layout(h,'force3')
%  view(3)
%%
A=D;
D=-1*(D-1);
for i=1:size(D,1)
    D(i,i)=0;
end
tree = linkage(D,'single','euclidean');
u=unique(names);
for k=1:length(u)
    ind=find(ismember(names,u{k}));
    truth_vec(ind)=k;
end
cluster_vec=cluster(tree,size(unique(names),2))

fscore(cluster_vec,truth_vec)
fscore_similar2(names',cluster_vec,truth_vec)

outperm = optimalleaforder(tree,D);
polardendrogram(tree,0,'colorthreshold','default','reorder',outperm,'labels',names);
view(2);
tightfig;
set(gcf,'Position',[100 100 800 800]);
save genelink4.mat D outperm names
saveas(gcf,'polardendrogram.pdf')
%%
u=unique(names);
U=zeros(length(u));
for i=1:length(u)
    for ii=1:length(u)
        j=sfind(names,u(i));
        jj=sfind(names,u(ii));
        U(i,ii)=nanmean(nanmean(A(j,jj)));
    end
    U(i,i)=1;
end
schemaball( U,u,[0 1 0],[1 0 0]); tightfig; tightfig2;
%saveas(gcf,'spiderweb.pdf');


%% spider digram with only certain names
%genes_not_to_include={'gcrR','341', '694', '1924', 'tatD'};
genes_not_to_include={'gcrR','341', '694', '1924', 'UA159', 'tatD'};

idx=[];
for i=genes_not_to_include
    idx=[idx,sfind(u,i)'];
end
idx=unique(idx);
idx2=1:length(u);
idx2(idx)=[];
U2=U;
U2(U2<.10)=0;
U3=sin(U2*pi/2)
geneid=[1.1930 1.3880 1.4270 1.6230 1.7810 0.4620 0.6950 0.8350 1.1280 1.9140 1.6720 0.9490 1.9160 1.9165 1.9170 0.0610 0.0630 1.9970 1.4270 0.5500 1.6910 1.3980 0.4860 0.8310 0.2450 0.4840 0.4890 0.4880 0.7530 0.9210];
[~,idx3]=sort(geneid);
schemaball4( U3(idx2(idx3),idx2(idx3)),u(idx2(idx3)),[0 .5 1],[1 .5 0]); tightfig; tightfig2;





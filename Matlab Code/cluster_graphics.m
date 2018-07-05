clear;
load clusters
load data 
close all

u = unique(names);
vec=zeros(1,length(names));
for k=1:length(u)
    ind=find(strcmp(names,u{k}),1);
    vec(ind)=k;
end
a=double(jet);
ind=round(vec/length(u)*63+1);
col=a(ind,:);



%%
c=clusters;
c(isnan(c))=0;
A=zeros(size(c,2)-2);
algo=size(c,1);
% m=mean(c(:,1).*c(:,2));
% syms x
% power=double(solve(algo==m^-x,x));

for jjj=1:algo
    c3=c(jjj,3:end);
    c1=c(jjj,1);
    c2=c(jjj,2);    
    for j=1:size(c,2)-2
        for jj=1:size(c,2)-2
            A(j,jj)=(c1*c2)^sensitivity*double(isequal(c3(j),c3(jj)))+A(j,jj);
        end
        A(j,j)=0;
    end
end
A=A/nanmax(A(:));
D=A;
%%
figure;
j=linspace(0,1,20);
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
saveas(gcf,[filepath filename '_network_plot.' figure_file_type])
 
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
saveas(gcf,[filepath filename '_polardendrogram.' figure_file_type])
%%
u=unique(names);
U=zeros(length(u));
for i=1:length(u)
    for ii=1:length(u)
        j=find(strcmp(names,u(i)),1);
        jj=find(strcmp(names,u(ii)),1);
        U(i,ii)=nanmean(nanmean(A(j,jj)));
    end
    U(i,i)=1;
end
schemaball3( U,u,[0 0 1],[0 1 0]);
c=colorbar;
set(c, 'ylim', [0 1])
c.Label.String='Relative Association';
c.Label.FontSize = 12;
saveas(gcf,[filepath filename '_spiderweb.' figure_file_type]);


% %% spider digram with only certain names
% genes_not_to_include=genes_to_exlude;
% 
% idx=[];
% for i=genes_not_to_include
%     s=strfind(names,i);
%     s(cellfun(@isempty,s))={0};
%     idx=[idx,find(cell2mat(s)>0)];
% end
% l=1:length(names);
% l(idx)=[];
% idx2=unique(l);
% U2=U;
% U2(U2<.10)=0;
% U3=sin(U2*pi/2);
% schemaball( U3(idx2,idx2),u(idx2),[0 .5 1],[1 .5 0]); tightfig; tightfig2;





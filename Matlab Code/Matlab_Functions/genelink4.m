clear all;
load clusters.mat
close all
clear genes
strain=strain_names3;
strain=replace(strain,'BHI','');
strain=replace(strain,'CSP','');
strain=replace(strain,' ','');
strain=replace(strain,'SMu','Smu');
strain=replace(strain,'Smu.','');
strain=replace(strain,'ftsQ','ftsQ/divIB');
strain=replace(strain,'IGR.1924','gcrR');
strain=replace(strain,'51','purK');
strain=replace(strain,'61','comR');
strain=replace(strain,'247','sufC');
strain=replace(strain,'286','nlmT');
strain=replace(strain,'485','liaF');
strain=replace(strain,'550','divIB');
strain=replace(strain,'915c','queF');
strain=replace(strain,'916c','queE');
strain=replace(strain,'917c','queD');
strain=replace(strain,'919c','queC');
strain=replace(strain,'921','rcrR');
strain=replace(strain,'949','clpX');
strain=replace(strain,'1128','ciaH');
strain=replace(strain,'1398','irvR');
strain=replace(strain,'1427','cdaR');
strain=replace(strain,'1489','lacX');
strain=replace(strain,'1541','pulA');
strain=replace(strain,'1675','metB');
strain=replace(strain,'1691','dltA');
strain=replace(strain,'1914c','cipB');
strain=replace(strain,'1916','comD');
strain=replace(strain,'1917','comE');
strain=replace(strain,'341','tatD');
strain=replace(strain,'919c','queC');
strain=replace(strain,'831','ltaS');
strain=replace(strain,'245','mecA');
strain=replace(strain,'833','rgpI');
strain=replace(strain,'753','pspC');
strain=replace(strain,'1826','yfbQ');
strain=replace(strain,'489','pnpB');
strain=replace(strain,'488','ppiB');
strain=replace(strain,'486','liaS');
strain=replace(strain,'484','pknB');
strain=replace(strain,'comX2','comX');
strain=replace(strain,'comR2','comR');
strain=replace(strain,'clpP2','clpP');
strain=replace(strain,'WT','UA159');
strain=regexprep(strain,'\([^)]*\)','');



u = unique(strain);
vec=zeros(1,length(strain));
for k=1:length(u)
    ind=sfind(strain,u{k});
    vec(ind)=k;
end
a=double(jet);
ind=round(vec/length(u)*64);
col=a(ind,:);



%%
c=clusters;
A=zeros(size(c,2)-2);
algo=size(c,1);
m=mean(c(:,1).*c(:,2));
syms x
power=double(solve(algo==m^-x,x));
for jjj=1:algo
    c3=c(jjj,3:end);
    c1=c(jjj,1);
    c2=c(jjj,2);
     %if c1 >.38
         %if c2 > .38
            %if length(unique(c3))>length(unique(c(:,3:end)))/1.1
                for j=1:size(c,2)-2
                    for jj=1:size(c,2)-2
                        A(j,jj)=(c1*c2)^1*double(isequal(c3(j),c3(jj)))+A(j,jj);
                    end
                    A(j,j)=0;
                end
            %end
          %end
      %end
end
A=A/max(A(:));
D=A;
%%
figure;
j=linspace(prctile(A(A>.85),0),1,20);
col2=winter;
nodenames=[];
for i=1:length(strain)
    ii=num2str(i);
    nodenames{i}=strcat(strain{i},'_',ii);
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
    %   G3=graph(A,strain1);
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
u=unique(strain_names3);
for k=1:length(u)
    ind=find(ismember(strain_names3,u{k}));
    truth_vec(ind)=k;
end
cluster_vec=cluster(tree,size(unique(strain_names3),2))

fscore(cluster_vec,truth_vec)
fscore_similar2(strain_names3',cluster_vec,truth_vec)

outperm = optimalleaforder(tree,D);
polardendrogram(tree,0,'colorthreshold','default','reorder',outperm,'labels',strain);
view(2);
tightfig;
set(gcf,'Position',[100 100 800 800]);
save genelink4.mat D outperm strain
saveas(gcf,'polardendrogram.pdf')
%%
u=unique(strain);
U=zeros(length(u));
for i=1:length(u)
    for ii=1:length(u)
        j=sfind(strain,u(i));
        jj=sfind(strain,u(ii));
        U(i,ii)=nanmean(nanmean(A(j,jj)));
    end
    U(i,i)=1;
end
schemaball( U,u,[0 1 0],[1 0 0]); tightfig; tightfig2;
%saveas(gcf,'spiderweb.pdf');


%% spider digram with only certain strain
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





%% heat maps

'heatmaps'
clear
close all
load data filename variables_to_include pvalue

% condition={'Acid','BHI','CSP','H2O2','O2','Paraquat','thermal'};
% for k=condition
eval(['load ' filename '_tnorm.mat tnorm cdf']);
%stoping point work on for loop below
for j={'control','Control','wildtype','Wild Type','WT','wt','wildtype'}'
    str=strcat('\textcolor{blue}{',j,'}');
    tnorm(:,1)=replace(tnorm(:,1),j,str);
end
try
    control=[];
    control=sfind(tnorm(2:end,1),'\textcolor{blue}');
end
% if strcmp(variables_to_include,'auto')==1 | isempty(variables_to_include)==1
%     variables=find(mean(cdf)<pvalue)+1;
% else
%     variables=sfind(tnorm(1,:),variables_to_include);
% end
% mat=cell2mat(tnorm(2:end,variables));
variables=2:size(tnorm,2);
mat=cell2mat(tnorm(2:end,variables));
mat(isnan(mat))=0;
if isempty(control)==0
    zmat=(mat-mat(control,:))./nanstd(mat);
else
    zmat=(mat-nanmean(mat))./nanstd(mat);
end
col=flipud(redgreencmap(256));
mcdf=mean(cdf);
labels=strcat(tnorm(1,variables)',' (P_{value}= ',num2str(mean(cdf)'),')');  
cgo=clustergram(flip(zmat',1),'rowlabels',flip(labels),'columnlabels',tnorm(2:end,1),'colormap',col,'cluster',2);
%addXLabel(cgo,num2str(mean(cdf)'),'position',[-2,15],'FontSize',11)
addTitle(cgo,strcat(filename,' Heatmap'))
% h=plot(cgo);
% set( findobj(gcf,'type','axes'),'FontSize',11,'FontAngle','italic','FontWeight','bold','outerposition',[ -0.03    0.02    1.04    0.7618])
% set(h,'position',[200,200,1200,600])
% a=annotation('textbox',[.0,.19,.1,.55],'string',num2str(mean(cdf)'),'FontSize',11,'FitBoxToText','on','FontWeight','bold');
% pv=uicontrol('style','text');
% set(pv,'test','position',[0,0,.1,.2])
% set(pv,num2str(mean(cdf)),'location',[0,0,.1,.2])
% saveas(gcf,[filename '_heatmap.pdf'])

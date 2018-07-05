% create program for bioscreen
%% import data
clear
%delete *.mat
if exist('Settings.txt','file')==0
     run Settings.m
else
    eval(fileread('Settings.txt'))
end
files=0;
[data,time,Data]=deal([]);
ls
filelist=dir([savepath '*' filetype])


%import data
for j={filelist.name}
    files=1+files
    if cell2mat(strfind(j,'xls')) ~=0
        [~, ~, data] = xlsread(j{:},data_sheet);
    elseif  cell2mat(strfind(j,'txt')) ~=0
        data = dlmread(j{:});
    elseif cell2mat(strfind(j,'csv')) ~=0
        data = csvread(j{:});
    else
        error= 'unknown file type; please specify .txt, .csv or xls file types'
    end
    if files>1 & size(data(start_row:end,time_column_number),1)~=size(Data,1)
        error= 'files do not have the same amount of time points; try importing files individualy'
    else
%        time=[time,cell2mat(data(start_row:end,time_column_number))];
%       fix - get different timing to work
        time=[cell2mat(data(start_row:end,time_column_number))];

        if iscell(data)==1
            data(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),data)) = {''};
            data=cell2mat(data(start_row:end,:));
        end
        blank=mean(data(:,blank_column_number),2);
        data2=data(:,max([time_column_number,blank_column_number])+1:end);
        Data=[Data,data2-blank];
    end
end
names=textscan(fopen(names_file),'%s');
names=names{:};
%reformat data
nor=number_of_replicates;
od=zeros(size(Data,2)/nor,nor,size(Data,1));
for i=1:size(Data,2)/nor
    for ii=1:nor
        od(i,ii,:)=Data(:,i+ii-1);
    end
end
idx=find(contains({'secs','mins','hours','days'},time_unit));
conversion=[1/60,1,60,24*60];
time=time*conversion(idx);
filename=filelist.name
filename=filename(1:end-4);
save data
%% figures of replicates
if strcmp(figures_of_fits,'yes')==1
    'fits w figures'
    clear
    close all
    load data
    parfor j=1:size(od,1)
        dcolor(10,0);
        figure(j); hold on;
        set(j,'Position',[100, 0, 1920, 1080]/1.5)
        for jj=1:number_of_replicates
            try
            [f1,f2,f3,d1,d2,d3,ma,mi,md,m1,m2,i,ii,iii,iiii,iiiii,timeexp,timed,odexp,odd,odm,tlag]=deal([]);
            subplot(nor,5,[1,6,11]); hold on;
            odr=squeeze(od(j,jj,:));
            plot(time,odr); title([names(j) ' od']); hold on;
            legend('Sample 1','Sample 2','Sample 3','Location', 'southeast'); ylabel('od'); xlabel('minuets')
            
            [ma,mi,md,m1,m2,i,ii,iii,iiii,iiiii,timeexp,timed,odexp,odd,odm]=important_points(time,odr);
            subplot(nor,5,jj*5); hold on;
            plot(time(1:end-1),diff(odr));
            plot(time(1:end-2),diff(odr,2));
            title(['sample ' num2str(jj)]); ylabel('\Delta od'); xlabel('minuets')
            legend('dy/dt','dy^2/dt^2')
            
            subplot(nor,5,jj*5-3); hold on;
            %growth
            [f1,f2,f3,pv1]=fitod(timeexp,odexp,.01,'y',j);
            subplot(nor,5,jj*5-2); hold on;
            %death
            if length(timed)> 10
                [d1,d2,d3,pv2]=fitod(timed,odd,-.001,'y',j);
            end
            %timelag
            tlag=timelag2(time,odr,f3);
            
            %fits
            subplot(3,5,jj*5-1); hold on;
            scatter(time,odr,7,'g');
            plot(timeexp,expf(timeexp,[f1,f2,f3]),'r',timed,expf(timed,[d1,d2,d3]),'r');
            plot(time,m1*(time-tlag)+f3,'k:'); plot(time,f3*ones(size(time)),'k:');
            scatter(tlag,f3,'bo');
            scatter(timeexp,odexp, 'b.'); scatter(timed,odd, 'b.');
            set(gca,'ylim',[0,max(odr)*1.1]);
            title({[' timelag = ' num2str(tlag)],['\chi^2 p_{value}= ',num2str(pv1/2+pv2/2)]});
            ylabel('od'); %xlabel('minuets')
            end
        end
        tightfig;
        try; names(j)={num2str(names{j})}; end;
        saveas(j,[names{j} '.' fit_figure_file_type])
        close(j)
    end
else
    'Fit figures option not specified. Set to "yes" in setting file to produce fit figures'
end
%% fits wo figure
'unnomilized fits'
clear
load data

%for k=1:size(names,1)

nc=length(names)
[ma,m1,ea,eb,ec,da,db,dc,p_val1,p_val2,tlag,tlag2,tlag3,c,h0,h1,i4]=deal(zeros(nc,3));

for j=1:length(names)
    for jj=1:nor
        odr=squeeze(od(j,jj,:));
        if max(odr)>lower_bound && max(odr)<upper_bound %something grew
            %find key points
            [ma(j,jj),mi,md,m1(j,jj),m2,i,ii,iii,i4(j,jj),iiiii,timeexp,timed,odexp,odd,odm,tlag3(j,jj)]=important_points(time,odr);
            %growth phase
            [ea(j,jj),eb(j,jj),ec(j,jj)]=fitod(timeexp,odexp,.01);
            %death phase
            [da(j,jj),db(j,jj),dc(j,jj)]=fitod(timed,odd,-.001);
            %timelag2
            [tlag2(j,jj),h0(j,jj)]=timelag2(time,odr,ec(j,jj));
            %curvature
            c(j,jj)=curvature(time,odr);
            h1(j,jj)=odr(i4(j,jj));
        end
    end
end
%timelag 2 & 3
tlag1=log(ea)./eb;
save([filename '_fits.mat'])
%end
%% bootstrap
'bootstrap'
clear
load data
% for k=1:size(names,1)
%     condition(k)
%import_thermal_stress
bootstrap_samples
bs=bootstrap_samples;
number_of_cultures=length(names)
nc=number_of_cultures;

tic
[ma,m1,i4,ea,eb,ec,ed,eres,da,db,dc,dres,tlag1,tlag2,tlag3,c,h0,h1]=deal(zeros(nc,3,bs));
parfor j=1:nc
    for jj=1:nor
        for jjj=1:bs
            
            %randomize
            r=sort(randi([1,size(od,3)],size(od,3),1));
            timer=time(r);
            odr=squeeze(od(j,jj,r));
            
        if max(odr)>lower_bound && max(odr)<upper_bound %something grew
                % try
                %find key points
                [ma(j,jj,jjj),mi,md,m1(j,jj,jjj),m2,i,ii,iii,i4(j,jj,jjj),iiiii,timeexp,timed,odexp,odd,odm,tlag3(j,jj,jjj)]=important_points(timer,odr);
                %exp
                [ea(j,jj,jjj),eb(j,jj,jjj),ec(j,jj,jjj)]=fitod(timeexp,odexp,.01);
                %end
                %try
                %death
                [da(j,jj,jjj),db(j,jj,jjj),dc(j,jj,jjj)]=fitod(timed,odd,-.001);
                % end
                % try
                %timelag2
               [tlag2(j,jj,jjj),h0(j,jj,jjj)]=timelag2(timer,odr,ec(j,jj,jjj));
                %curvature
                c(j,jj,jjj)=curvature(timer,odr);
                h1(j,jj,jjj)=odr(i4(j,jj,jjj));
                % end
            end
        end
    end
end
tlag1=log(ea)./eb;
toc
save([filename '_bootstrap.mat'])
%eval(['save ',condition{k},'_bootstrap.mat']);
%end
%% create nomilized table via bootstap fits
'nomilized table'
load data
tic
bs=bootstrap_samples
if bs>0
    eval(['load ' filename '_bootstrap.mat ma m1 i4 ea eb ec da db dc tlag1 tlag2 tlag3 c h0 h1 names']);
    tbs=cat(4, ma, m1, i4, ea, eb, ec, da, db, dc, tlag1, tlag2, tlag3, h0, h1, c);
    tbs(abs(tbs)>1e5)=NaN;
    tbs(tbs==0)=NaN;
    stdv=squeeze(nanstd(tbs,1,3));
    stdv(stdv==0)=1e-5;
end
clearvars ma m1 ea eb ec da db dc tlag1 tlag2 tlag3 h0 h1 c
file=[filename '_fits.mat'];
eval(['load ' filename '_fits.mat ma m1 i4 ea eb ec da db dc tlag1 tlag2 tlag3 c h0 h1']);
tfit=cat(3, ma, m1, i4, ea, eb, ec, da, db, dc, tlag1, tlag2, tlag3, h0, h1, c);
tbs(abs(tbs)>1e5)=NaN;
tfit(abs(tfit)>1e5)=NaN; tfit(tbs==0)=NaN;
headers={'Strain' 'Carrying capacity \newline(~ #cells)' 'Maxmium \newlineGrowth Rate \newline(1/min)' 'Maxmium \newlineGrowth Delay \newline(min)' 'A \newline(Growth prefactor) \newline~#cells' 'Average \newlineGrowth rate \newline(1/min)' 'Growth constant \newline~ #cells' 'Death prefactor \newline~ #cells' 'Death rate \newline(1/min)' 'Death constant \newline~# cells' 'Lag time 1 \newline(min)' 'Lag time 2 \newline(min)' 'Lag time 3 \newline(min)' 'curvature' 'inflection point' 'od at tlag2'}
tnorm=[];
for j=1:length(names)
    %idx=sfind(names,u(j),'exact');
    tfit_i=tfit(j,:,:);
    if bs>0
        varbs=stdv(j,:,:).^2;
    end
    varfit=((tfit_i-nanmean(tfit_i))./nanstd(tfit_i)).^2;
    if strcmpi(filter_outliers,'no')==1 && bs==0
        tnorm(j,:)=nanmean(tfit_i,2);
    elseif strcmpi(filter_outliers,'no')==1 && bs>0
        numerator=nansum(tfit_i./varfit,2);
        denominator=nansum(1./varfit,2);
        tnorm(j,:)=numerator./denominator;
    elseif strcmpi(filter_outliers,'yes')==1 && bs==0
        numerator=nansum( tfit_i./varbs,2 );
        denominator=nansum( 1./varbs,2 );
        tnorm(j,:)=numerator./denominator;
    elseif strcmpi(filter_outliers,'yes')==1 && bs>0
        numerator=nansum( tfit_i./varbs./varfit,2 );
        denominator=nansum( 1./varbs./varfit,2 );
        tnorm(j,:)=numerator./denominator;
    else
        'unkown normilization options: set "filter_outliers=" to yes or no and set bootstrap_samples to a number grater than zero'
    end
    %tnorm(j,jj)=nansum(tfit(j,:,jjj)./stdv(j,:,jjj).^2)/nansum(1./stdv(j,:,jjj).^2);
end
tnorm=num2cell(tnorm);
tnorm=[names,tnorm];
tnorm=[headers;tnorm];

%%pvlaues
for i=1:size(tfit,1)
    for ii=1:size(tfit,2)
        for iii=1:size(tfit,3)
            chi_square(i,ii,iii)=(tfit(i,ii,iii)-cell2mat(tnorm(i+1,iii+1))).^2./nanstd(tfit(:,ii,iii)).^2;
        end
    end
end
chi_square3=squeeze(nansum(chi_square,2));
cdf=chi2cdf(chi_square3,3);
save([filename '_tnorm.mat']);
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
if strcmp(variables_to_include,'auto')==1 | isempty(variables_to_include)==1
    variables=find(mean(cdf)<pvalue)+1;
else
    variables=sfind(tnorm(1,:),variables_to_include);
end
mat=cell2mat(tnorm(2:end,variables));
mat(isnan(mat))=0;
if isempty(control)==0
    zmat=(mat-mat(control,:))./nanstd(mat);
else
    zmat=(mat-nanmean(mat))./nanstd(mat);
end
col=flipud(redgreencmap(256));
%zmat(isnan(zmat))=0;
h=plot(clustergram(flip(zmat',1),'rowlabels',flip(tnorm(1,variables)),'columnlabels',tnorm(2:end,1),'colormap',col,'cluster',2));
ax = findobj(gcf,'type','axes');
set(ax,'FontSize',11,'FontAngle','italic','FontWeight','bold','outerposition',[ -0.03    0.02    1.04    0.7618])
saveas(gcf,[filename '_heatmap.pdf'])
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
        [~,~,~,~,~,~,~,f,~,fsimilar,T]=dendqual5(t3,z,metrics{jj},methods{j},r,names,[filename,'_cluster_vecs.csv']);
        r=1+r
        dlmwrite([filename '_cluster_measures.csv'],[f,fsimilar,T'],'delimiter', ',' ,'-append')
        clusters(r,:)=[f,fsimilar,T'];
        %end
    end
end
save clusters.mat clusters names t3 f fsimilar filename
toc
'average fscore='
mean(clusters(:,1:2))
%% cluster graphics
run cluster_graphics
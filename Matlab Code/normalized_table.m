%% create nomilized table via bootstap fits
clear
'nomilized table'
load data
tic
bs=bootstrap_samples
if bs>0
    eval(['load ' filename '_bootstrap.mat ma m1 i1 ea eb ec da db dc tlag1 tlag2 tlag3 c h0 h1 names']);
    tbs=cat(4, ma, m1, i1, ea, eb, ec, da, db, dc, tlag1, -tlag2, tlag3, h0, h1, c);
    tbs(abs(tbs)>1e5)=NaN;
    tbs(tbs==0)=NaN;
    stdv=squeeze(nanstd(tbs,1,3));
    stdv(stdv==0)=1e-5;
end
clearvars ma m1 ea eb ec da db dc tlag1 tlag2 tlag3 h0 h1 c
file=[filelist.name(1,:) '_fits.mat'];
eval(['load ' filename '_fits.mat ma m1 i1 ea eb ec da db dc tlag1 tlag2 tlag3 c h0 h1']);
tfit=cat(3, ma, m1, i1, ea, eb, ec, da, db, dc, tlag1, -tlag2, tlag3, h0, h1, c);
tfit(abs(tfit)>1e5)=NaN; 
headers={'Strain' 'Carrying capacity (~ #cells)' 'Maxmium Growth Rate (1/min)' 'Maxmium Growth Delay (min)' 'A (Growth prefactor) ~#cells' 'Average Growth rate (1/min)' 'Growth constant ~ #cells' 'Death prefactor ~ #cells' 'Death rate (1/min)' 'Death constant ~# cells' 'Lag time 1 (min)' 'Lag time 2 (min)' 'Lag time 3 (min)' 'Curvature' 'OD @ timelag' 'OD @ max Growth'}
tnorm=[];
for j=1:length(names)
    %idx=sfind(names,u(j),'exact');
    tfit_i=tfit(j,:,:);
    if bs>0
        tbs(abs(tbs)>1e5)=NaN;
        tfit(tbs==0)=NaN;
        varbs=stdv(j,:,:).^2;
    end
    varfit=((tfit_i-nanmean(tfit_i))./nanstd(tfit_i)).^2;
    if strcmpi(filter_outliers,'no')==1 && bs==0
        tnorm(j,:)=nanmean(tfit_i,2);
    elseif strcmpi(filter_outliers,'no')==1 && bs>0
        numerator=nansum( tfit_i./varbs./2 );
        denominator=nansum( 1./varbs./2 );
        tnorm(j,:)=numerator./denominator;
    elseif strcmpi(filter_outliers,'yes')==1 && bs==0
        numerator=nansum( tfit_i./varfit,2 );
        denominator=nansum( 1./varfit,2 );
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

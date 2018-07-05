function [m,mm,mmm,m1,m2,i,ii,iii,iiii,iiiii,timeexp,timed,odexp,odd,odm,tlag]=important_points(timer,odr)


% odr=squeeze(od(j,jj,:));
% timer=time(j)';
[m,i]=max(odr); %carrying capacity
[mm,ii]=min(odr); %carrying capacity offset
[mmm,iii]=min(abs(odr-(m/2-mm/2))); %mid point between mid and max

odm=odr(iii:i);
[m2,iiiii]=min(diff(odm,2));
[m1,iiii]=max(diff(odr)./diff(timer));
%[~,iiiiii]=min(diff(odr,2));


odexp=odr(ii:iiii);
timeexp=timer(ii:iiii);
timed=[];
odd=[];
if length(odr)-i-7>1
    odd=odr(round((length(odr)-i)/2)+i:end);
    timed=timer(round((length(odr)-i)/2)+i:end);
% else
%     odd=odr(end);
%     timed=timer(end);
end
[ut,tidx]=unique(timeexp);
if length(ut)>2
    tt=timeexp(1):timeexp(end);
    eexp=ppval(spline(ut,odexp(tidx)),tt);
    [~,tidx2]=min((eexp-min(eexp)-.1).^2);
    tlag=-tt(tidx2);
else
    tlag=NaN;
end
end
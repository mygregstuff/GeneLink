function [timelag,ho]=timelag2(time,od,offset)
[time,i]=unique(time);
od=od(i);
[m,ii]=max(diff(od)./diff(time));
%y=mx+b
b=od(ii)-m*time(ii);
x=[min(time):.1:max(time)];
res=@(x)(m*x(1)+b-offset).^2;
%opts=optimset('display','off');
timelag=fminsearch(res,200);%,opts);
od_sp=spline(time,od,0:1:time(end));
ho=od_sp(round(timelag));

% opts =  optimset('display','off');
% initalguess=-222;
% j=0;
% while diff(time(i-j:i+1))<15
%     j=j+1;
% end      
% timelag = fminsearch('res2',initalguess,opts,time(i-j:i+1),od(i-j:i+1),m,offset);
 end

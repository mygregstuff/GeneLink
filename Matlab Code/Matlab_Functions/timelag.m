function delay_time=timelag(time,od,offset)
[utime,idx]=unique(time);
sp=spline(utime,od(idx),0:time(end));
[~,delay_time]=min((sp-offset-1).^2);
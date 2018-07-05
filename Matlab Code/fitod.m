function [fa,fb,fc,p_value]=fitod(t,y,u,err_est,drawplot)

% options=optimset('MaxFunEvals',1e5, 'MaxIter',5e4, 'Tolfun', 1e-4, 'TolX',1e-4);
if size(t)==size(y) & length(t)>4
    t1 = t(1) ; t2 = t(end);
    y1 = y(1) ; y2 = y(end);
    a = (y1-y2)/(exp(u*t1)-exp(u*t2));
    c = y1-a*exp(u*t1);
    f1 = [a u c];
    
    opts =  optimset('display','off');
    [fitpar,sumsq] = fminsearch('res3',f1,opts,t,y);
    
    if exist('err_est','var')==1
        chi_square=sumsq/var(y);
        p_value=chi2cdf(chi_square,length(y));
    end
    
    % fitpar=fminsearch('res3',f1,options,t,y);
    %fa = abs(fitpar(1));
    fa = fitpar(1);
    fb = fitpar(2);
    %fc = abs(fitpar(3));
    fc = fitpar(3);
    
    % resnorm=res(fitpar(4),t,y);
    
    if exist('drawplot','var')==1
        figure(drawplot);
        yguess=expf(t,f1);
        y_fit=expf(t,[fa,fb,fc]);
        semilogy(t,y,'b.',t,yguess,'b',t,y_fit,'r');
        title3=strcat('a=', num2str(fa),  '  b=', num2str(fb), '  c=', num2str(fc));
        title2=strcat('a=', num2str(a),  '  b=', num2str(u), '  c=', num2str(c));
        title({title2;title3}); ylabel('od'); xlabel('minuets');
        ylabel('od'); xlabel('minuets');
        legend({'Initial guees','Fit'})
    end    
else
%     "error: check the size of the data"
%     t
%     y
[fa,fb,fc,p_value]=deal(NaN);
end
end

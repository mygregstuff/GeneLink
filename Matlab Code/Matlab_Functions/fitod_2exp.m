function [fa,fb,sumsq]=fitod_2exp(t,y,u,sigma,drawplot)

% options=optimset('MaxFunEvals',1e5, 'MaxIter',5e4, 'Tolfun', 1e-4, 'TolX',1e-4);

t1 = t(1) ; t2 = t(end);
y1 = y(1) ; y2 = y(end);
 a = (y1-y2)/(exp(u*t1)-exp(u*t2));
 c = y1-a*exp(u*t1); 
f1 = [u c];

opts =  optimset('display','off'); 
[fitpar, sumsq] = fminsearch('res2',f1,opts,t,y);
 
fa = fitpar(1);
fb = fitpar(2);


if exist('sigma','var')==1
chi_square=sumsq/var(y);
p_val=chi2cdf(chi_square,length(y));
end

if exist('drawplot','var')==1
    figure(drawplot);
    yguess=expf(t,f1);
    y_fit=expf(t,[exp(-1*fa*fb),fb,fc]);
    semilogy(t,y,'b.',t,yguess,'b',t,y_fit,'r');
    title3=strcat('a=', num2str(fa),  '  b=', num2str(fb), '  c=', num2str(fc));
    title2=strcat('a=', num2str(a),  '  b=', num2str(u), '  c=', num2str(c));
    title({title2;title3}); ylabel('od'); xlabel('minuets');
    ylabel('od'); xlabel('minuets');
end
end

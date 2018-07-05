function sumsq = res2(tlag,t,y,m,offset)%,drawplot)

% function sumsq = explfit(A,t,y)
%  model y = A1 exp (- A2 t) + A3
% 
%  where A = abs(A)
y_fit = m*(t+tlag)+offset;
sumsq = sum(((y_fit-y)).^2);


% if exist('drawplot')==1
%     figure(drawplot);
%     semilogy(t,y,'.',t,y_fit,'r');
% end
end
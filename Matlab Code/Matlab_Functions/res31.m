function sumsq = res3(A,t,y)%,drawplot)

% function sumsq = explfit(A,t,y)
%  model y = A1 exp (- A2 t) + A3
% 
%  where A = abs(A)

% A(1)=abs(A(1));
% A(3)=abs(A(3));

%y_fit = A(1) * exp(A(2)*t) + A(3);
y_fit = exp(A(2)*(t-A(1))) + A(3);

sumsq = sum(((y_fit-y)./y).^2);


% if exist('drawplot')==1
%     figure(drawplot);
%     semilogy(t,y,'.',t,y_fit,'r');
% end
end
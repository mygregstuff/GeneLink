function [idx,dc,d]=tvsm(X,varargin)
% idx=find(strcmp(varargin,'metric1'));
% metric1=varargin(idx+1);
% idx=find(strcmp(varargin,'metric2'));
% metric2=varargin(idx+1);
% idx=find(strcmp(varargin,'method'));
% method=varargin(idx+1);
% pd=pdist(X,metric1{:});
% tree=linkage(X,method{:},metric2{:});
% idx=find(strcmp(varargin,'Criteria'));
% Criteria=varargin(idx+1);
% idx=find(strcmp(varargin,'Transformation'));
% Transformation=varargin(idx+1);
pd=pdist(X);
tree=linkage(X,'weighted');
idx=optimalleaforder(tree,pd);
dc=sum(sqrt(sum(diff(X(idx,:)).^2,2)));


if exist('varargin','var')==1
    if strcmp(varargin,'figure')==1
        %scatter(X(:,1),X(:,2)); hold on
        plot(X(idx,1),X(idx,2),'r','LineWidth',1);
        axis equal
        grid on
        legend(['D = ',num2str(dc)])
        title 'Clustering Solution'
    elseif sum(strcmp(varargin,'true'))==1  %% find the true solution
        close gcf
        npoints=size(X,1);
        p=gpu(perms(1:npoints));
        X=gpu(X);
        x=X(:,1);
        y=X(:,2);
        dx=diff(x(p)');
        dy=diff(y(p)');
        D=sum(sqrt(dx.^2+dy.^2));
        %sequntail points from the permiations
        [d,idx]=min(D);
        idx=p(idx,:);
        if dc-d>.01 && sum(strcmp(varargin,'figure'))==1
            figure('pos',[100,100,1200,600]);
            subplot(1,2,1)
            scatter(X(:,1),X(:,2)); hold on
            axis equal
            grid on
            plot(X(idx_c,1),X(idx_c,2),'r','LineWidth',1);
            legend(['D = ',num2str(d_c)])
            title 'Clustering Solution'
            
            subplot(1,2,2)
            scatter(X(idx,1),X(idx,2),'ro'); hold on
            axis equal
            grid on
            plot(x(idx),y(idx),'b-');
            legend(['D = ',num2str(d)])
            title 'True Solution'
            tightfig
        end
    end
end
end
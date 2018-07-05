function [idx,d,d_c]=tvsm(X,varargin)
X=gpu(X); %move array over to gpu
tree=gpu(linkage(X))
pd=gpu(pdist(X))
idx_c=optimalleaforder(tree,pd,'Criteria','group')
d_c=sum(sqrt(sum(diff(X(idx_c,:)).^2,2)))
idx=idx_c;
if exist('varargin','var')==1
    if strcmp(varargin,'figure')==1
        scatter(X(:,1),X(:,2)); hold on
        axis equal
        grid on
        plot(X(idx_c,1),X(idx_c,2),'r','LineWidth',1);
        legend(['D = ',num2str(d_c)])
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
        if d_c-d>.01 && sum(strcmp(varargin,'figure'))==1
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
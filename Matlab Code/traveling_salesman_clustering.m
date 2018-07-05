%traveling salesman proplem
%steps
%determine clusters
%connect n-1 nearest negbors in clusters (n is cluster size)
%connect clusters endpoints stariung with nearst clusters and working up
%check to make sure theres no overlaps
clear
close all
profile on
npoints=10;
X=rand(npoints,2);
% d=zeros(npoints);
% x=X(:,1);
% y=X(:,2);
% scatter(x,y, 'bo'); hold on;
% axis equal
% grid on
% for i=1:npoints
%     for j=1:npoints
%         dx=x(j)-x(i);
%         dy=y(j)-y(i);
%         d(i,j)=sqrt(dx^2+dy^2); %distance matrix
%     end
% end
%% radial intial guess
% scatter(mean(x),mean(y),'go'); hold on;%find the centroid
% for i=1:npoints
%     radi(i)=norm(x(i)-mean(x),y(i)-mean(y));%find the radi
%     ang(i)=atan((y(i)-mean(y))/(x(i)-mean(x)));  %find the radial angle
% end
% [~,idx]=sort(ang);
% plot(x(idx),y(idx),'g:'); hold on;
%connect points based on radial angle nad radius
%% shortest paths guessing

% [~,idx]=sort(d(:));
% for i=1:npoints
%     idx2=idx(i+npoints);
%     j(i)=mod(idx2,npoints);
%     k(i)=(idx2-j(i))/npoints+1;
%     if j(i)==0
%         j(i)=npoints;
%     end
%     %if sum(j(i)==j)<2 || sum(k(i)==k)<2
%         plot([x(j(i)),x(k(i))],[y(j(i)),y(k(i))],'g-')
%     %end
% end
%% baysian 
X=gpu(X);
x=X(:,1);
y=X(:,2);
c=nchoosek(1:npoints,3);
for i=1:size(c,1)
    p=perms(c(i,:));
    p=p(1:3,:); %for speed
    dx=diff(x(p)');
    dy=diff(y(p)');
    d=sum(sqrt(dx.^2+dy.^2));
    [~,idx]=min(d);
    bestpaths(i,:)=p(idx,:); 
end
paths=nchoosek(1:npoints,2);
freq=zeros(nchoosek(npoints,2),1);
for i=1:size(paths,1)
    p1=(bestpaths==paths(i,1));
    p2=(bestpaths==paths(i,2));
    s1=sum(p1(:,1).*p2(:,2));
    s2=sum(p1(:,2).*p2(:,1));
    s3=sum(p1(:,2).*p2(:,3));
    s4=sum(p1(:,3).*p2(:,2));
    freq(i)=freq(i)+s1+s2+s3+s4;
end
idx=find(freq==max(freq));
p=paths(idx,:);
length(p)
dx=diff(x(p)');
dy=diff(y(p)');
d_peicewise=sum(sqrt(dx.^2+dy.^2))

%% clutering
tree=linkage([x,y],'average');
pd=pdist(X);
idx=optimalleaforder(tree,pd);
d_cluster=sum(sqrt(sum(diff(X(idx,:)).^2,2)))

%% find the true solution
npoints=size(X,1);
p=gpu(perms(1:npoints));
X=gpu(X);
x=X(:,1);
y=X(:,2);
dx=diff(x(p)');
dy=diff(y(p)');
d=sum(sqrt(dx.^2+dy.^2));
%sequntail points from the permiations
d_true=min(d)
profile viewer

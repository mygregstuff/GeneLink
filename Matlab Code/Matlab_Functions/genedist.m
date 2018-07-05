function [gd,v]=genedist(Dsquare,vec)
[gd,v]=deal(0);
for i=unique(vec)
    ind=find(vec==i);
    lind=length(ind);
    if lind>1
        Md=Dsquare(ind,ind)./mean(Dsquare(:));
        gd=sum(Md(:))/lind/(lind-1)+gd;
        listd=diff(ind);
        v=poissfit(listd)+v;
    end
end
end
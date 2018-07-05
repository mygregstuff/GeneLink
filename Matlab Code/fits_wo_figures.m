%% fits wo figure
'unnomilized fits'
clear
load data

%for k=1:size(names,1)

nc=length(names);
[ma,m1,ea,eb,ec,da,db,dc,p_val1,p_val2,tlag,tlag2,tlag3,c,h0,h1,i1]=deal(zeros(nc,3));

parfor j=1:length(names)
    for jj=1:nor
        odr=squeeze(od(j,jj,:));
        if max(odr)>lower_bound && max(odr)<upper_bound %something grew
            %find key points
            [ma(j,jj),mi,md,m1(j,jj),m2,i,ii,iii,i1(j,jj),iiiii,timeexp,timed,odexp,odd,odm,tlag3(j,jj)]=important_points(time,odr);
            %growth phase
            [ea(j,jj),eb(j,jj),ec(j,jj)]=fitod(timeexp,odexp,.01);
            %death phase
            [da(j,jj),db(j,jj),dc(j,jj)]=fitod(timed,odd,-.001);
            %timelag2
            [tlag2(j,jj),h0(j,jj)]=timelag2(time,odr,ec(j,jj));
            %curvature
            c(j,jj)=curvature(time,odr);
            h1(j,jj)=odr(i1(j,jj));
        end
    end
end
%timelag 2 & 3
tlag1=log(ea)./eb;
save([filename '_fits.mat'])
%end

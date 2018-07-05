%% bootstrap
'bootstrap'
clear
load data
% for k=1:size(names,1)
%     condition(k)
%import_thermal_stress
bootstrap_samples
bs=bootstrap_samples;
number_of_cultures=length(names)
nc=number_of_cultures;

tic
[ma,m1,i1,ea,eb,ec,ed,eres,da,db,dc,dres,tlag1,tlag2,tlag3,c,h0,h1]=deal(zeros(nc,3,bs));
parfor j=1:nc
    for jj=1:nor
        for jjj=1:bs
            
            %randomize
            r=sort(randi([1,size(od,3)],size(od,3),1));
            timer=time(r);
            odr=squeeze(od(j,jj,r));
            
        if max(odr)>lower_bound && max(odr)<upper_bound %something grew
                try
                %find key points
                [ma(j,jj,jjj),mi,md,m1(j,jj,jjj),m2,i,ii,iii,i1(j,jj,jjj),iiiii,timeexp,timed,odexp,odd,odm,tlag3(j,jj,jjj)]=important_points(timer,odr);
                end
                try
                %exp
                [ea(j,jj,jjj),eb(j,jj,jjj),ec(j,jj,jjj)]=fitod(timeexp,odexp,.01);
                end
                try
                %death
                [da(j,jj,jjj),db(j,jj,jjj),dc(j,jj,jjj)]=fitod(timed,odd,-.001);
                end
                try
                %timelag2
               [tlag2(j,jj,jjj),h0(j,jj,jjj)]=timelag2(timer,odr,ec(j,jj,jjj));
                %curvature
                end
                try
                c(j,jj,jjj)=curvature(timer,odr);
                h1(j,jj,jjj)=odr(i1(j,jj,jjj));
                end
            end
        end
    end
end
tlag1=log(ea)./eb;
toc
save([filename '_bootstrap.mat'])
%eval(['save ',condition{k},'_bootstrap.mat']);
%end

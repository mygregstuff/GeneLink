%% figures of replicates
if strcmp(figures_of_fits,'yes')==1
    'fits w figures'
    clear
    close all
    load data
    parfor j=1:size(od,1)
        time2=time;
        dcolor(10,0);
        figure(j); hold on;
        set(j,'Position',get( groot,'Screensize')*nor/3)
        for jj=1:number_of_replicates
            [f1,f2,f3,d1,d2,d3,ma,mi,md,m1,m2,i,ii,iii,iiii,iiiii,timeexp,timed,odexp,odd,odm,tlag]=deal([]);
            subplot(nor,5,1:5:nor*5-4); hold on;
            odr=squeeze(od(j,jj,:));
            plot(time,odr); title([names(j) ' od']); hold on;
            legend('Sample 1','Sample 2','Sample 3','Location', 'southeast'); ylabel('od'); xlabel('minuets')
            %try
            [ma,mi,md,m1,m2,i,ii,iii,iiii,iiiii,timeexp,timed,odexp,odd,odm]=important_points(time,odr);
            subplot(nor,5,jj*5); hold on;
            plot(time2(1:end-1),diff(odr));
            plot(time2(2:end-1),diff(odr,2));
            title(['sample ' num2str(jj)]); ylabel('\Delta od'); xlabel('minuets')
            legend('dy/dt','dy^2/dt^2')
            % end
            subplot(nor,5,jj*5-3); hold on;
            %growth
            [f1,f2,f3,pv1]=fitod(timeexp,odexp,.01,'y',j);
            subplot(nor,5,jj*5-2); hold on;
            %death
            if length(timed)> 10
                [d1,d2,d3,pv2]=fitod(timed,odd,-.001,'y',j);
            end
            %timelag
            tlag=timelag2(time,odr,f3);
            
            %fits
            subplot(nor,5,jj*5-1); hold on;
            scatter(time,odr,7,'g');
            plot(timeexp,expf(timeexp,[f1,f2,f3]),'r',timed,expf(timed,[d1,d2,d3]),'r');
            plot(time,m1*(time-tlag)+f3,'k:'); plot(time,f3*ones(size(time)),'k:');
            scatter(tlag,f3,'bo');
            scatter(timeexp,odexp, 'b.'); scatter(timed,odd, 'b.');
            set(gca,'ylim',[0,max(odr)*1.1]);
            title({[' timelag = ' num2str(tlag)],['\chi^2 p_{value}= ',num2str(pv1/2+pv2/2)]});
            ylabel('od'); %xlabel('minuets')
            %end
        end
        tightfig;
        try; names(j)={num2str(names{j})}; end;
        saveas(j,[names{j} '.' fit_figure_file_type])
        close(j)
    end
else
    'Fit figures option not specified. Set to "yes" in setting file to produce figures'
end

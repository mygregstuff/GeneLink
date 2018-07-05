function [d]=similargenedist(Dsquare,genes)

genegroup={{'Smu.550 (ftsQ); +','Smu.485 (liaF); +','Smu.1427c'};...
    {'Smu.1916 (comD); +','comE','Smu.1916 (comD); -','IGR61 (comS?)','Smu.1128 (ciaH); -','Smu.61 (comR); -','Smu.286 (nlmT); -','IGR61','IGR1914 (cipB promoter)','IGR1914 (cipB promoter)','Smu.1914c (cipB); -'};...
    {'Smu.1909c','Smu.1909c; -','Smu.1910c'};...
    {'Smu.914c; +','Smu.915c (queF); -','Smu.916c','Smu.916c (queE); +','Smu.916c (queE); -','Smu.919c (queC); -','Smu.917c'}};

d=0;
md=mean(Dsquare(:));
for i=1:length(genegroup)
    c=find(ismember(genes,{genegroup{i}{:}}))';
    l=length({genegroup{i}{:}});
    for j=c
        for jj=c
            if ismember( genes(j), genes(jj))==1
                d=Dsquare(j,jj)*l/length(c)/md+d;
            end
        end
    end
%     d=sum(D2(:))/mean(Dsquare(:))/(numel(D2)-sum(c))+d;
end
end
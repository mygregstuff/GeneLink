function f=fscore_similar2(genes,cluster_vec,truth_vec)
norm = zeros(size(unique(cluster_vec)));
[TP,TN,FP,FN] = deal(zeros(length(unique(cluster_vec)),length(unique(truth_vec))));
% genegroup={{'Smu.550 (ftsQ); +','Smu.485 (liaF); +','Smu.1427c'};...
%     {'Smu.1916 (comD); +','comE','Smu.1916 (comD); -','IGR61 (comS?)','Smu.1128 (ciaH); -','Smu.61 (comR); -','Smu.286 (nlmT); -','IGR61','IGR1914 (cipB promoter)','IGR1914 (cipB promoter)','Smu.1914c (cipB); -','rcrR promoter'};...
%     {'Smu.1909c','Smu.1909c; -','Smu.1910c'};...
%     {'Smu.914c; +','Smu.915c (queF); -','Smu.916c','Smu.916c (queE); +','Smu.916c (queE); -','Smu.919c (queC); -','Smu.917c'}};
genegroup={{'Smu.550','Smu.485','Smu.1427'};...
    {'comD','IGR61','Smu.1128','Smu.61','Smu.286','IGR61','IGR1914','Smu.1914','rcrR promoter'};...
    {'Smu.1909','Smu.1910'};...
    {'Smu.914','Smu.915','Smu.916','Smu.919','Smu.917'}
    {'Smu.63','Smu.1946','Smu.928','Smu.1945'}};
for i=genegroup'
    subind=[];
    for ii=i{:}
        subind=[subind,find(~cellfun(@isempty,strfind(genes,ii)))];        
    end   
    %subind=find(ismember(genes,i{:}));
    if isempty(subind)==0    
        cluster_vec(ismember(cluster_vec,subind))=subind(1);
        truth_vec(ismember(truth_vec,subind))=subind(1);
    end
end
for i=unique(cluster_vec')
    norm(i) = sum(cluster_vec==i)/length(cluster_vec);
    for ii=unique(truth_vec)
        TP(i,ii) = sum(truth_vec(cluster_vec==i)==ii);
        FP(i,ii) = sum(truth_vec(cluster_vec==i)~=ii);
        FN(i,ii) = sum(truth_vec(cluster_vec~=i)==ii);
        TN(i,ii) = sum(truth_vec(cluster_vec~=i)~=ii);
    end
end
f = nansum(max( 2*TP ./ (2*TP + FP + FN) ,[],2) .* norm);
% [i,ii]=meshgrid(unique(cluster_vec'),unique(truth_vec));
% TP(i,ii)=sum(truth_vec(cluster_vec==i)==ii);
% TP=bsxfun(@(i,ii)sum(truth_vec(cluster_vec==i)==ii),unique(cluster_vec'),unique(truth_vec))
end
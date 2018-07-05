function f=fscore(cluster_vec,truth_vec)
norm = zeros(size(unique(cluster_vec)));
[TP,TN,FP,FN] = deal(zeros(length(unique(cluster_vec)),length(unique(truth_vec))));
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
% ind=bsxfun(@find,[unique(cluster_vec)]==cluster_vec);
% ind=bsxfun(@find,truth_vec(ind)==unique(truth_vec));

% TP=bsxfun(@iseqaule,ind,truth_vec)
%  TP(:,:)=arrayfun(@(i,ii)sum(truth_vec(cluster_vec==i)==ii),1:4,1:4)
% TP(:,:)=arrayfun(@(ii)sum(truth_vec(ismember(ind,unique(cluster_vec')))==ii),unique(truth_vec))

end
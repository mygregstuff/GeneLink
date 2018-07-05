function [idx,numbers]=ssort(cell)
n=[regexp(cell,'\d*','Match')];
[n,idx] = sort([n{:}]);
numbers=cell2mat(n');
end
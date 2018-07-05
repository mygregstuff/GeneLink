function [vector,idx]=str2num2(cell_of_strings_with_numbers)
n=regexp(cell_of_strings_with_numbers,'\d*','Match');
n=str2double([n{:}]);
[vector,idx]=sort(n);
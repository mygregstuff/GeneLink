function idx=sfinde(cell_array,string
a=strfind(cell_array,string);
a(cellfun(@isempty,a))={0};
a=cell2mat(a);
idx=find(a==1);
end

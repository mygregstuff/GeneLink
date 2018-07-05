%% import data
clear
delete *.mat
run Settings.m
files=0;
[Data,Time]=deal([]);
filelist=dir(filetype)
for j={filelist.name}
    files=1+files
    if cell2mat(strfind(j,'xls')) ~=0
        [~, ~, data] = xlsread(j{:},data_sheet);
    elseif  cell2mat(strfind(j,'txt')) ~=0
        data = dlmread(j{:});
    elseif cell2mat(strfind(j,'csv')) ~=0
        data = csvread(j{:});
    else
        error= 'unknown file type; please specify .txt, .csv or xls file types'
    end
    if files>1 & size(data,1)~=Data
        error= 'files do not have the same amount of time points; try importing files individualy'
    else
        time_units={'Time','time','Minuests','minuets','Hours','hours','Seconds','seconds'}; jj=1;
        while isempty(find(strcmp(data,time_units(jj)),1))==1  
            jj+jj+1;
        end
        [tidx1,tidx2]=find(strcmp(data,time_units(jj)),1);
        Time=cell2mat(data(tidx1+1,tidx2));
        control={'Blank','blank','Water','water','Control','control'}; jj=1;
        while isempty(find(strcmp(data,control(jj)),1))==1  
            jj+jj+1;
        end        
        [bidx1,bidx2]=find(strcmp(data,control(jj)),1);
        blank=mean(cell2mat(data(bidx1+1,bidx2)));
        data2=data(tidx1+1:end,max([bidx2,tidx2])+1:end);
        data2(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),data2)) = {''};
        data2=cell2mat(data2);
        Data=[Data,data2-blank];
    end
end
names=textscan(fopen(names_file),'%s');
names=names{:};
%reformat data
nor=number_of_replicates;
od=zeros(size(Data,2)/nor,nor,size(Data,1));
for i=1:size(Data,2)/nor
    for ii=1:nor
        od(i,ii,:)=Data(:,i+ii-1);
    end
end
time=Time*60*60;
% stopping point - fix fielname to now include extension
i=cell2mat(strfind(cellstr(filelist.name),'.'));
filename=filelist.name(1:i-1)
save data

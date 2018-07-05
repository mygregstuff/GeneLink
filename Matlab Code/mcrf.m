function []=mcrf(script_to_be_compiled)
tic

string={"mcc -o script_to_be_compiled -W main:script_to_be_compiled -T link:exe -d './script_to_be_compiled/for_testing' -v 'script_to_be_compiled.m' -a '../*.m'"}
newstring=strrep(string{:},'script_to_be_compiled',script_to_be_compiled)
eval(newstring{:})

string={"copyfile ../*.txt script_to_be_compiled/for_testing/"}
newstring=strrep(string{:},'script_to_be_compiled',script_to_be_compiled)
eval(newstring{:})

string={"copyfile ../*.xls script_to_be_compiled/for_testing/"}
newstring=strrep(string{:},'script_to_be_compiled',script_to_be_compiled)
eval(newstring{:})

string={"system('sh ./script_to_be_compiled/for_testing/run_script_to_be_compiled.sh /Applications/MATLAB/MATLAB_Runtime/v93')"}
newstring=strrep(string{:},'script_to_be_compiled',script_to_be_compiled)
eval(newstring{:})


toc
end

% percent signs comment out sections of code and is not used 
% settings
data_select='manual'
filetype='.xls' %specify file type
start_row=2 %row which the data starts
bootstrap_samples=200
figures_of_fits='no'
lower_bound=.3 %discards data with an optical density maximum below a given value
%helps remove curves where nothing grew
upper_bound=1.5 %discards data with an optical density maximum above a given value
%helps remove curves with optical aberrations nothing grew
figure_file_type='pdf' %png, svg, pdf, jpeg... etc
filter_outliers='yes'
different_conditions={'controle','contorle+codition'}
time_unit='days' % sec, min, hours, or days.
number_of_replicates=3 % users can either use one file or have a list of names of for each separate
strains_to_exclude={} %list of names from names.txt to not in include in analysis... If you know some strains are bad
variables_to_include='auto'%must be atlaeast one of following:'auto' or {'carrying capacity' 'Maximum' 'Lag time' 'Growth prefactor' 'Growth rate' 'Growth constant' 'Death prefactor' 'Death rate' 'Death constant'}
%example: variables_to_include={'carrying capacity' 'Maximum' 'Lag time' 'Growth rate' 'Death prefactor'}
pvalue=.25 %dont include growth parameters larger than a set pvalue, part of with variables_to_include='auto'
number_of_independent_varibles=1
%sensitivity of adjacency matrix
sensitivity=.5;
%time = 1 independent variable; 
%time and 2 concentrations of glucose = 2 independent variables
%time and 3 concentrations of glucose = 3 independent variables
%f_simialr:file of similar genes used for calculating the fsimialr scores
%simply just change it
Quit_after_run='no' % yes automatically closes the program once default figures are produced

try; nor=round(length(names)/length(unique(names))); end;%number of replicates
try; nor=round(length(names)/length(unique([names{:}]))); end;%number of replicates
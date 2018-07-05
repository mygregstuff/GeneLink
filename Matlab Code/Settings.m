% percent signs comment out sections of code and is not used 
% settings
filetype='*.xls' %specify file type
data_sheet=''
start_row=2 %row which the data starts
time_column_number=1 %colum number of the time
blank_column_number=2
bootstrap_samples=20
figures_of_fits='no'
lower_bound=.3 %discards data with an optical desity maxmium below a given value
%helps remove curves where nothing grew
upper_bound=1.5 %discards data with an optical desity maxmium above a given value
%helps remove curves with optical aberations nothing grew
fit_figure_file_type='pdf' %png, svg, pdf, jpeg... etc
filter_outliers='yes'
different_conditions={'bhi','bhi+csp'}
time_unit='days' % sec, min, hours, or days.
number_of_replicates=3 % users can either use one file or have a list of names of for each sperate
names_file='names.txt'; % strain names should contain replicates if present. % Wild type, wild type wild type , non wild type, etc....
strains_to_exclude={} %list of names from names.txt to not in include in analysis... If you know some strains are bad
variables_to_include='auto'%must be atlaeast one of following:'auto' or {'carrying capacity' 'Maxmium' 'Lag time' 'Growth prefactor' 'Growth rate' 'Growth constant' 'Death prefactor' 'Death rate' 'Death constant'}
%example: variables_to_include={'carrying capacity' 'Maxmium' 'Lag time' 'Growth rate' 'Death prefactor'}
pvalue=.25 %dont include growth paramaters larger than a set pvalue, part of with variables_to_include='auto'
number_of_independent_varibles=1
%sensitivty of adjacency matrix
sensitivity=.5;
%time = 1 independent variable; 
%time and 2 concentrations of gluecose = 2 independent varibles
%time and 3 concentrations of gluecose = 3 independent varibles
%f_simialr:file of similar genes used for calulating the fsimialr scores
%simply just change it
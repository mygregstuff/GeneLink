% percent signs comment out sections of code and is not used 
% settings
filetype='*.xls'
time_column_number=1
blank_column_number=2
start_row=2
data_sheet=''
bootstrap_samples=3
figures_of_fits='no'
%png, svg, pdf, jpeg... etc
fit_figure_file_type='pdf'
filter_outliers='yes'
different_conditions={'bhi','bhi+csp'}
number_of_replicates=3
% uses can eith use one file or have a list of names of for each sperate
% file. Names should contain replicates if replacest are present 
% Wild type, wild type wild type , non wild type, etc....
names_file='names.txt';
%must be one of the strian names in the file names.txt
strains_to_exclude={}
%must be atlaeast one of following:'auto' or {'carrying capacity' 'Maxmium' 'Lag time' 'Growth prefactor' 'Growth rate' 'Growth constant' 'Death prefactor' 'Death rate' 'Death constant'}
%variables_to_include={'carrying capacity' 'Maxmium' 'Lag time' 'Growth rate' 'Death prefactor'}
variables_to_include='auto'
pvalue=.25
%time = 1 independent variable; 
%time and 2 concentrations of gluecose = 2 independent varibles
%time and 3 concentrations of gluecose = 3 independent varibles
number_of_independent_varibles=1
%f_simialr:file of similar genes used for calulating the fsimialr scores
%simply just change it
%Produce aggregated ied table with variables of interest, and perform
%Poisson generalized linear mixed effects model.
%
%Different sections of this script produce statistical analysis for the 
%effect of flicker on interictal epileptiform discharge counts.
%
%Script is designed to run in sections; see the comment header for purpose
%of each script section (generate and save a simplified/summed data table,
%run statistical analysis, or save statistical output as a data structure 
%for later figure editing)
%2024/02/27

%% Set dataset directory and dependent repos:
[root_dir,repo]=define_flicker_root_dir;
for field_name=fieldnames(repo)'
    addpath(genpath(repo.(field_name{:})));
end

%note: can use function fetch_flicker_subjectIDs to get list of sessions
%and their details.

%% Define paths:

output_dir=[root_dir '/stg-analyses/ied-analysis'];
source_data=[root_dir '/stg-analyses/NatureComm2024-figures/Source Data']; %where source data for figures is stored

%% Produce table aggregating ied counts and variables of interest:
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end
produce_flicker_ied_tbl(root_dir,output_dir);

%% Perform GLM analyses and save output:

ied_GLM;


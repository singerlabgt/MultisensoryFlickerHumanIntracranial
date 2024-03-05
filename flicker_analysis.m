%This script analyses data from the flickerneuro or flickerfreq task from begining to end.
%2024/02/25

%% Set dataset directory and dependent repos:
[root_dir,repo]=define_flicker_root_dir;
for field_name=fieldnames(repo)'
    addpath(genpath(repo.(field_name{:})));
end

%note: can use function fetch_flicker_subjectIDs to get list of sessions
%and their details.

%% manually define parameters (including which sections of the script to run):

% params.subjectID={''}; %subject ID
% params.task={''}; %task (flickerneuro or flickerfreq)
% params.ses={''}; %session number
% params.fetchTrials=0;
% params.IDnoisych=0;
% params.ied_chExclude=0;
% params.preproc=0;
% params.runPSD=0;
% params.eval_degree_ent=0;
% params.eval_plv=0;
% params.plot_ent_by_stimfreq=0;
% params.echo_analysis=0;
% params.plot_ent_3d=0;
% params.single_unit_analysis=0;

%% set additional parameters and paths:

%the following sections are not run if task is flickerneuro
if strcmp(params.task,'flickerneuro')
    params.eval_plv=0;
    params.plot_ent_by_stimfreq=0;
elseif strcmp(params.task,'flickerfreq')
    params.single_unit_analysis=0;
end

%set paths:
fnames=struct();
fnames.subjectID=params.subjectID{:}; %subject ID
fnames.task=params.task{:}; %task
fnames.ses=params.ses{:}; %session number
fnames.root_dir=define_flicker_root_dir; %directory where data is
fnames.rawdata_folder=[fnames.root_dir '/stg-raw/sub-' fnames.subjectID '/task-' fnames.task '/ses-' fnames.ses]; %not included in publication
fnames.LFP_filename=[fnames.rawdata_folder '/sub-' fnames.subjectID '_stg-raw_task-' fnames.task '_ses-' fnames.ses '_nat-ieeg.EDF']; %not included in publication
fnames.exclude_ch=[fnames.root_dir '/stg-raw/sub-' fnames.subjectID '/anat/sub-' fnames.subjectID '_stg-raw_desc-chExclude_nat-list.txt']; %not inlucded in publication
fnames.preprocdata_folder=[fnames.root_dir '/stg-preproc/sub-' fnames.subjectID '/task-' fnames.task '/ses-' fnames.ses]; %preprocessed data folder
fnames.noisy_ch=[fnames.preprocdata_folder '/sub-' fnames.subjectID '_stg-preproc_task-' fnames.task '_ses-' fnames.ses '_desc-noisy-ch_nat-list.csv']; %where list of noisy channels is stored
fnames.trials_filename=[fnames.preprocdata_folder '/sub-' fnames.subjectID '_stg-preproc_task-' fnames.task '_ses-' fnames.ses '_nat-beh.mat']; %structure with information about trials
fnames.preproc_LFPdata=[fnames.preprocdata_folder '/LFP/static_ent/']; %preprocessed LFP data
fnames.analysis_folder=[fnames.root_dir '/stg-analyses/task-' fnames.task '/sub-' fnames.subjectID '/ses-' fnames.ses]; %analyzed data folder

%% run PSD and save output:
if params.runPSD
    static_ent_PSDanalysis(fnames);
end

%% calculate significance and amplitude of modulation for all electrodes:
if params.eval_degree_ent
    evaluate_ent_degree_relpower(fnames);
end

%% calculate phase-locking value significance and amplitude between flicker and LFP:
if params.eval_plv
    evaluate_ssep_plv(fnames);
end

%% run analysis on persistence of oscillatory activity:
if params.echo_analysis
    echo_analysis(fnames);
end

%% draw degree of entrainment in function of stimulation frequency, for each channel:
if params.plot_ent_by_stimfreq
    plot_mod_by_stimfreq(fnames,'Laplacian');
end

%% (optional) draw individual's modulation on 3D model of subject's brain, and save output:
if params.plot_ent_3d
    plot_ind_entrainment_3d(fnames,'Laplacian');
end

%% if single units, show modulation of single units
%Note: need to have separately preprocessed single units before-hand.
if params.single_unit_analysis
    if exist([fnames.preprocdata_folder '/single-unit-preproc'],'dir')
        unit=importdata([fnames.preprocdata_folder '/single-unit-preproc/sub-' fnames.subjectID '_stg-preproc_task-' fnames.task '_ses-' fnames.ses '_nat-unit-data.mat'],'unit'); %load single unit data
        if ~isempty(unit) %if there are units, analyze their modulation
            modulation_singleunits(fnames);
        end
    end
end


%This script analyzes data from the spep task from begining to end.
%2024/02/25

%% Set dataset directory and dependent repos:
[root_dir,repo]=define_flicker_root_dir;
for field_name=fieldnames(repo)'
    addpath(genpath(repo.(field_name{:})));
end

%note: can use function fetch_flicker_subjectIDs to get list of sessions
%and their details.

%% manually define parameters (including which sections of the script to run):

% params.subjectID={''};
% params.ses={''};
% params.task={'spep'};
% params.fetchTrials=0;
% params.IDnoisych=0;
% params.preproc=0;
% params.runERP=0;
% params.eval_degree_spep=0;
% params.plot_spep_group_3d=0;

%% set some paths:

fnames=struct();
fnames.subjectID=params.subjectID{:};
fnames.task=params.task{:};
fnames.ses=params.ses{:};
fnames.root_dir=define_flicker_root_dir;
fnames.rawdata_folder=[fnames.root_dir '/stg-raw/sub-' fnames.subjectID '/task-spep/ses-' fnames.ses];
fnames.LFP_filename=[fnames.rawdata_folder '/sub-' fnames.subjectID '_stg-raw_task-spep_ses-' fnames.ses '_nat-ieeg.EDF'];
fnames.exclude_ch=[fnames.root_dir '/stg-raw/sub-' fnames.subjectID '/anat/sub-' fnames.subjectID '_stg-raw_desc-chExclude_nat-list.txt'];
fnames.preprocdata_folder=[fnames.root_dir '/stg-preproc/sub-' fnames.subjectID '/task-spep/ses-' fnames.ses];
fnames.noisy_ch=[fnames.preprocdata_folder '/sub-' fnames.subjectID '_stg-preproc_task-spep_ses-' fnames.ses '_desc-noisy-ch_nat-list.csv'];
fnames.trials_filename=[fnames.preprocdata_folder '/sub-' fnames.subjectID '_stg-preproc_task-spep_ses-' fnames.ses '_nat-beh.mat'];
fnames.preproc_LFPdata=[fnames.preprocdata_folder '/LFP/spep_erp/'];
fnames.analysis_folder=[fnames.root_dir '/stg-analyses/task-spep/sub-' fnames.subjectID '/ses-' fnames.ses];

%% run ERP analysis and save output:
if params.runERP
    spep_ERPanalysis(fnames);
end

%% establish which contacts have SPEP response, and amplitude and timing of response:
if params.eval_degree_spep
    evaluate_SPEP_degree(fnames);
end

%% plot SPEP significance and amplitude on 3D plot:
if params.plot_spep_group_3d
    plot_group_spep_3d(fnames,{fnames.subjectID},'Laplacian');
end

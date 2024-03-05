%Simulate 40Hz flicker data based on superposition hypothesis, and run PSD analysis.
%2024/02/25

%% Set dataset directory and dependent repos:
[root_dir,repo]=define_flicker_root_dir;
for field_name=fieldnames(repo)'
    addpath(genpath(repo.(field_name{:})));
end

%note: can use function fetch_flicker_subjectIDs to get list of sessions
%and their details.

%% manually define subjectID:
% params.subjectID={''};

%% set paths:

fnames=struct();
fnames.subjectID=params.subjectID{:};
fnames.task='spep';
fnames.ses='01';
fnames.root_dir=define_flicker_root_dir;
fnames.rawdata_folder=[fnames.root_dir '/stg-raw/sub-' fnames.subjectID '/task-spep/ses-' fnames.ses];
fnames.LFP_filename=[fnames.rawdata_folder '/sub-' fnames.subjectID '_stg-raw_task-spep_ses-' fnames.ses '_nat-ieeg.EDF'];
fnames.exclude_ch=[fnames.root_dir '/stg-raw/sub-' fnames.subjectID '/anat/sub-' fnames.subjectID '_stg-raw_desc-chExclude_nat-list.txt'];
fnames.preprocdata_folder=[fnames.root_dir '/stg-preproc/sub-' fnames.subjectID '/task-spep/ses-' fnames.ses];
fnames.noisy_ch=[fnames.preprocdata_folder '/sub-' fnames.subjectID '_stg-preproc_task-spep_ses-' fnames.ses '_desc-noisy-ch_nat-list.csv'];
fnames.trials_filename=[fnames.preprocdata_folder '/sub-' fnames.subjectID '_stg-preproc_task-spep_ses-' fnames.ses '_nat-beh.mat'];
fnames.preproc_LFPdata=[fnames.preprocdata_folder '/LFP/spep_erp/'];
fnames.analysis_folder=[fnames.root_dir '/stg-analyses/model-superposition/sub-' fnames.subjectID '/ses-' fnames.ses];

%% produce model data:
superposition_modeling_data(fnames);

%% process modeled data like we process flicker data:
load([fnames.analysis_folder '/sub-' fnames.subjectID '_stg-analysis_model-superposition_ses-01_nat-trialsdata-refLaplacian.mat'],'trials','data');
cfg=[];
cfg.demean         = 'yes'; %apply baseline correction %SHOULD WE DO BASELINE CORRECTION?
cfg.baselinewindow = 'all'; %in seconds- the default ('all') is the complete trial %DO WE NEED TO CHANGE THIS?
cfg.lpfilter       = 'yes'; %lowpass filter
cfg.lpfreq         = 300; %lowpass frequency- same as Abby
cfg.hpfilter       = 'yes'; %highpass filter %FIGURE OUT WHY HIGH PASS FILTER DOESN'T WORK
cfg.hpfreq         = 2; %highpass frequency
data_ref_preproc=ft_preprocessing(cfg,data);

%run PSD on processed model data:
%requires to use MATLAB's dpss (not FieldTrip's)
PSD_results_ref_preproc=run_PSD_flicker(data_ref_preproc,trials);

%create output folder if does not exist:
if ~exist([fnames.analysis_folder '/modeldata_entrainment-PSDs'],'dir')
    mkdir([fnames.analysis_folder '/modeldata_entrainment-PSDs']);
end
plot_depthelectrode_PSDs(PSD_results_ref_preproc,[fnames.analysis_folder '/modeldata_entrainment-PSDs']);

%save output:
save([fnames.analysis_folder '/sub-' fnames.subjectID '_stg-analysis_model-superposition_ses-' fnames.ses '_nat-psd-refLaplacian.mat'],'PSD_results_ref_preproc');

%% evaluate significance and amplitude of modulation for simulated data:
temp=load([fnames.analysis_folder '/sub-' fnames.subjectID '_stg-analysis_model-superposition_ses-' fnames.ses '_nat-psd-refLaplacian.mat'],'PSD_results_ref_preproc');
PSD_results=temp.PSD_results_ref_preproc; %choose which version  of PSD results you want to use

conditions_of_interest=["40Hz-V","40Hz-AV","40Hz-A"];
control_condition='Baseline';
pvalue_table=zeros(length(PSD_results.label),length(conditions_of_interest));
zscore_table=zeros(length(PSD_results.label),length(conditions_of_interest));
for i=1:size(zscore_table,2) %for each condition
    freq_interest=str2num(regexprep(conditions_of_interest{i},'Hz.+',''));
    freq_interest_index=find(PSD_results.data{1,strcmp(PSD_results.condition,conditions_of_interest{i})}{3}==freq_interest);
    
    for ch=1:size(zscore_table,1) %for each recording channel
        stim_values=[];
        baseline_values=[];
        for tr=1:size(PSD_results.data{ch,strcmp(PSD_results.condition,conditions_of_interest{i})}{1},1) %for however many number of trials of given condition there are
            current_stim_value=PSD_results.data{ch,strcmp(PSD_results.condition,conditions_of_interest{i})}{1}(tr,freq_interest_index);
            current_baseline_value=PSD_results.data{ch,strcmp(PSD_results.condition,control_condition)}{1}(tr,freq_interest_index);
            
            stim_values=[stim_values current_stim_value];
            baseline_values=[baseline_values current_baseline_value];
        end
        zscore_table(ch,i)=(mean(stim_values)/mean(baseline_values))-1;
        pvalue_table(ch,i)=pval_randomshuffle([stim_values' baseline_values'],10000);
    end
end

%save output:
zscore_table=array2table(zscore_table,'RowNames',PSD_results.label,'VariableNames',conditions_of_interest); %make matrix into table
writetable(zscore_table,[fnames.analysis_folder,'/LFP_zscore_table_ref-Laplacian.csv'],'WriteRowNames',1);

pvalue_table=array2table(pvalue_table,'RowNames',PSD_results.label,'VariableNames',conditions_of_interest); %make matrix into table
writetable(pvalue_table,[fnames.analysis_folder,'/LFP_pvalue_table_ref-Laplacian.csv'],'WriteRowNames',1);

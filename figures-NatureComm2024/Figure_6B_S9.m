%ied_gen_figs.m
%
%generate figures made in lou_stats_v3 using saved estimates and p-values
%in data structure IED_results.mat
%2024/02/26

%% load data and define outputs:

root_dir=define_flicker_root_dir;
results_folder=[root_dir '/stg-analyses/ied-analysis'];
ied_results=importdata([results_folder '/IED_struct.mat']);
ied_by_sess=importdata([results_folder '/IED_by_session.mat']);

source_data=[root_dir '/stg-analyses/NatureComm2024-figures/Source Data']; %where source data for figures is stored
figures_dir=[root_dir '/stg-analyses/NatureComm2024-figures/'];

%% aggregate data of interest:

results={};

%effects of any flicker on all patients:
results(end+1,:)={'all',ied_results.soz_type.pct_est(strcmp(ied_results.soz_type.labels,'All')),ied_results.soz_type.pct_ci(strcmp(ied_results.soz_type.labels,'All'),:),ied_results.soz_type.pvals(strcmp(ied_results.soz_type.labels,'All')),ied_results.soz_type.df(strcmp(ied_results.soz_type.labels,'All'))};
results(end+1,:)=repmat({''},1,5);

%interaction of flicker modality and region of interest:
results(end+1,:)={'visual-flicker-in-visual',ied_results.region_modality.pct_est(strcmp(ied_results.region_modality.labels,'Visual'),strcmp(ied_results.region_modality.modalities,'Visual')),[ied_results.region_modality.pct_ci_lower(strcmp(ied_results.region_modality.labels,'Visual'),strcmp(ied_results.region_modality.modalities,'Visual')) ied_results.region_modality.pct_ci_upper(strcmp(ied_results.region_modality.labels,'Visual'),strcmp(ied_results.region_modality.modalities,'Visual'))],ied_results.region_modality.pvals(strcmp(ied_results.region_modality.labels,'Visual'),strcmp(ied_results.region_modality.modalities,'Visual')),ied_results.region_modality.df(strcmp(ied_results.region_modality.labels,'Visual'),strcmp(ied_results.region_modality.modalities,'Visual'))};
results(end+1,:)={'visual-flicker-in-auditory',ied_results.region_modality.pct_est(strcmp(ied_results.region_modality.labels,'Audio'),strcmp(ied_results.region_modality.modalities,'Visual')),[ied_results.region_modality.pct_ci_lower(strcmp(ied_results.region_modality.labels,'Audio'),strcmp(ied_results.region_modality.modalities,'Visual')) ied_results.region_modality.pct_ci_upper(strcmp(ied_results.region_modality.labels,'Audio'),strcmp(ied_results.region_modality.modalities,'Visual'))],ied_results.region_modality.pvals(strcmp(ied_results.region_modality.labels,'Audio'),strcmp(ied_results.region_modality.modalities,'Visual')),ied_results.region_modality.df(strcmp(ied_results.region_modality.labels,'Audio'),strcmp(ied_results.region_modality.modalities,'Visual'))};
results(end+1,:)={'visual-flicker-in-MTL',ied_results.region_modality.pct_est(strcmp(ied_results.region_modality.labels,'MTL'),strcmp(ied_results.region_modality.modalities,'Visual')),[ied_results.region_modality.pct_ci_lower(strcmp(ied_results.region_modality.labels,'MTL'),strcmp(ied_results.region_modality.modalities,'Visual')) ied_results.region_modality.pct_ci_upper(strcmp(ied_results.region_modality.labels,'MTL'),strcmp(ied_results.region_modality.modalities,'Visual'))],ied_results.region_modality.pvals(strcmp(ied_results.region_modality.labels,'MTL'),strcmp(ied_results.region_modality.modalities,'Visual')),ied_results.region_modality.df(strcmp(ied_results.region_modality.labels,'MTL'),strcmp(ied_results.region_modality.modalities,'Visual'))};
results(end+1,:)={'visual-flicker-in-PFC',ied_results.region_modality.pct_est(strcmp(ied_results.region_modality.labels,'PFC'),strcmp(ied_results.region_modality.modalities,'Visual')),[ied_results.region_modality.pct_ci_lower(strcmp(ied_results.region_modality.labels,'PFC'),strcmp(ied_results.region_modality.modalities,'Visual')) ied_results.region_modality.pct_ci_upper(strcmp(ied_results.region_modality.labels,'PFC'),strcmp(ied_results.region_modality.modalities,'Visual'))],ied_results.region_modality.pvals(strcmp(ied_results.region_modality.labels,'PFC'),strcmp(ied_results.region_modality.modalities,'Visual')),ied_results.region_modality.df(strcmp(ied_results.region_modality.labels,'PFC'),strcmp(ied_results.region_modality.modalities,'Visual'))};
results(end+1,:)=repmat({''},1,5);
results(end+1,:)={'audio-flicker-in-visual',ied_results.region_modality.pct_est(strcmp(ied_results.region_modality.labels,'Visual'),strcmp(ied_results.region_modality.modalities,'Audio')),[ied_results.region_modality.pct_ci_lower(strcmp(ied_results.region_modality.labels,'Visual'),strcmp(ied_results.region_modality.modalities,'Audio')) ied_results.region_modality.pct_ci_upper(strcmp(ied_results.region_modality.labels,'Visual'),strcmp(ied_results.region_modality.modalities,'Audio'))],ied_results.region_modality.pvals(strcmp(ied_results.region_modality.labels,'Visual'),strcmp(ied_results.region_modality.modalities,'Audio')),ied_results.region_modality.df(strcmp(ied_results.region_modality.labels,'Visual'),strcmp(ied_results.region_modality.modalities,'Audio'))};
results(end+1,:)={'audio-flicker-in-auditory',ied_results.region_modality.pct_est(strcmp(ied_results.region_modality.labels,'Audio'),strcmp(ied_results.region_modality.modalities,'Audio')),[ied_results.region_modality.pct_ci_lower(strcmp(ied_results.region_modality.labels,'Audio'),strcmp(ied_results.region_modality.modalities,'Audio')) ied_results.region_modality.pct_ci_upper(strcmp(ied_results.region_modality.labels,'Audio'),strcmp(ied_results.region_modality.modalities,'Audio'))],ied_results.region_modality.pvals(strcmp(ied_results.region_modality.labels,'Audio'),strcmp(ied_results.region_modality.modalities,'Audio')),ied_results.region_modality.df(strcmp(ied_results.region_modality.labels,'Audio'),strcmp(ied_results.region_modality.modalities,'Audio'))};
results(end+1,:)={'audio-flicker-in-MTL',ied_results.region_modality.pct_est(strcmp(ied_results.region_modality.labels,'MTL'),strcmp(ied_results.region_modality.modalities,'Audio')),[ied_results.region_modality.pct_ci_lower(strcmp(ied_results.region_modality.labels,'MTL'),strcmp(ied_results.region_modality.modalities,'Audio')) ied_results.region_modality.pct_ci_upper(strcmp(ied_results.region_modality.labels,'MTL'),strcmp(ied_results.region_modality.modalities,'Audio'))],ied_results.region_modality.pvals(strcmp(ied_results.region_modality.labels,'MTL'),strcmp(ied_results.region_modality.modalities,'Audio')),ied_results.region_modality.df(strcmp(ied_results.region_modality.labels,'MTL'),strcmp(ied_results.region_modality.modalities,'Audio'))};
results(end+1,:)={'audio-flicker-in-PFC',ied_results.region_modality.pct_est(strcmp(ied_results.region_modality.labels,'PFC'),strcmp(ied_results.region_modality.modalities,'Audio')),[ied_results.region_modality.pct_ci_lower(strcmp(ied_results.region_modality.labels,'PFC'),strcmp(ied_results.region_modality.modalities,'Audio')) ied_results.region_modality.pct_ci_upper(strcmp(ied_results.region_modality.labels,'PFC'),strcmp(ied_results.region_modality.modalities,'Audio'))],ied_results.region_modality.pvals(strcmp(ied_results.region_modality.labels,'PFC'),strcmp(ied_results.region_modality.modalities,'Audio')),ied_results.region_modality.df(strcmp(ied_results.region_modality.labels,'PFC'),strcmp(ied_results.region_modality.modalities,'Audio'))};
results(end+1,:)=repmat({''},1,5);
results(end+1,:)={'audiovisual-flicker-in-visual',ied_results.region_modality.pct_est(strcmp(ied_results.region_modality.labels,'Visual'),strcmp(ied_results.region_modality.modalities,'AV')),[ied_results.region_modality.pct_ci_lower(strcmp(ied_results.region_modality.labels,'Visual'),strcmp(ied_results.region_modality.modalities,'AV')) ied_results.region_modality.pct_ci_upper(strcmp(ied_results.region_modality.labels,'Visual'),strcmp(ied_results.region_modality.modalities,'AV'))],ied_results.region_modality.pvals(strcmp(ied_results.region_modality.labels,'Visual'),strcmp(ied_results.region_modality.modalities,'AV')),ied_results.region_modality.df(strcmp(ied_results.region_modality.labels,'Visual'),strcmp(ied_results.region_modality.modalities,'AV'))};
results(end+1,:)={'audiovisual-flicker-in-auditory',ied_results.region_modality.pct_est(strcmp(ied_results.region_modality.labels,'Audio'),strcmp(ied_results.region_modality.modalities,'AV')),[ied_results.region_modality.pct_ci_lower(strcmp(ied_results.region_modality.labels,'Audio'),strcmp(ied_results.region_modality.modalities,'AV')) ied_results.region_modality.pct_ci_upper(strcmp(ied_results.region_modality.labels,'Audio'),strcmp(ied_results.region_modality.modalities,'AV'))],ied_results.region_modality.pvals(strcmp(ied_results.region_modality.labels,'Audio'),strcmp(ied_results.region_modality.modalities,'AV')),ied_results.region_modality.df(strcmp(ied_results.region_modality.labels,'Audio'),strcmp(ied_results.region_modality.modalities,'AV'))};
results(end+1,:)={'audiovisual-flicker-in-MTL',ied_results.region_modality.pct_est(strcmp(ied_results.region_modality.labels,'MTL'),strcmp(ied_results.region_modality.modalities,'AV')),[ied_results.region_modality.pct_ci_lower(strcmp(ied_results.region_modality.labels,'MTL'),strcmp(ied_results.region_modality.modalities,'AV')) ied_results.region_modality.pct_ci_upper(strcmp(ied_results.region_modality.labels,'MTL'),strcmp(ied_results.region_modality.modalities,'AV'))],ied_results.region_modality.pvals(strcmp(ied_results.region_modality.labels,'MTL'),strcmp(ied_results.region_modality.modalities,'AV')),ied_results.region_modality.df(strcmp(ied_results.region_modality.labels,'MTL'),strcmp(ied_results.region_modality.modalities,'AV'))};
results(end+1,:)={'audiovisual-flicker-in-PFC',ied_results.region_modality.pct_est(strcmp(ied_results.region_modality.labels,'PFC'),strcmp(ied_results.region_modality.modalities,'AV')),[ied_results.region_modality.pct_ci_lower(strcmp(ied_results.region_modality.labels,'PFC'),strcmp(ied_results.region_modality.modalities,'AV')) ied_results.region_modality.pct_ci_upper(strcmp(ied_results.region_modality.labels,'PFC'),strcmp(ied_results.region_modality.modalities,'AV'))],ied_results.region_modality.pvals(strcmp(ied_results.region_modality.labels,'PFC'),strcmp(ied_results.region_modality.modalities,'AV')),ied_results.region_modality.df(strcmp(ied_results.region_modality.labels,'PFC'),strcmp(ied_results.region_modality.modalities,'AV'))};
results(end+1,:)=repmat({''},1,5);

%modulated vs non-modulated channels:
results(end+1,:)={'non-mod-chs',ied_results.flicker_strength.pct_est(1),ied_results.flicker_strength.pct_ci(:,1)',ied_results.flicker_strength.pvals(1),ied_results.flicker_strength.df(1)}; %ADD DEGREES OF FREEDOM
results(end+1,:)={'mod-chs',ied_results.flicker_strength.pct_est(2),ied_results.flicker_strength.pct_ci(:,2)',ied_results.flicker_strength.pvals(2),ied_results.flicker_strength.df(2)}; %ADD DEGREES OF FREEDOM
results(end+1,:)=repmat({''},1,5);

%experiment 1 flicker effects on IEDs by condition, on whole brain:
flickerneuro_conditions={'5.5Hz-V','5.5Hz-AV','5.5Hz-A','40Hz-V','40Hz-AV','40Hz-A','80Hz-V','80Hz-AV','80Hz-A','R-V','R-AV','R-A'};
for cond=flickerneuro_conditions
    results(end+1,:)={['exp1_' cond{:} '_whole-brain'],ied_results.condition_exp1.pct_est(strcmp(ied_results.condition_exp1.labels,cond)),ied_results.condition_exp1.pct_ci(strcmp(ied_results.condition_exp1.labels,cond),:),ied_results.condition_exp1.pvals(strcmp(ied_results.condition_exp1.labels,cond)),ied_results.condition_exp1.df(strcmp(ied_results.condition_exp1.labels,cond))};
end
results(end+1,:)=repmat({''},1,5);

%experiment 2 flicker effects on IEDs by condition, on whole brain:
flickerfreq_conditions={'5.5Hz','8Hz','11Hz','14Hz','17Hz','20Hz','23Hz','26Hz','29Hz','32Hz','35Hz','38Hz','40Hz','42Hz','45Hz','48Hz','51Hz','54Hz','57Hz','63Hz','66Hz','69Hz','72Hz','75Hz','78Hz','80Hz','R'};
for cond=strcat(flickerfreq_conditions,'-V')
    results(end+1,:)={['exp2_' cond{:} '_whole-brain'],ied_results.condition_exp2.pct_est(strcmp(ied_results.condition_exp2.labels,cond)),ied_results.condition_exp2.pct_ci(strcmp(ied_results.condition_exp2.labels,cond),:),ied_results.condition_exp2.pvals(strcmp(ied_results.condition_exp2.labels,cond)),ied_results.condition_exp2.df(strcmp(ied_results.condition_exp2.labels,cond))};
end
results(end+1,:)=repmat({''},1,5);
for cond=strcat(flickerfreq_conditions,'-A')
    results(end+1,:)={['exp2_' cond{:} '_whole-brain'],ied_results.condition_exp2.pct_est(strcmp(ied_results.condition_exp2.labels,cond)),ied_results.condition_exp2.pct_ci(strcmp(ied_results.condition_exp2.labels,cond),:),ied_results.condition_exp2.pvals(strcmp(ied_results.condition_exp2.labels,cond)),ied_results.condition_exp2.df(strcmp(ied_results.condition_exp2.labels,cond))};
end
results(end+1,:)=repmat({''},1,5);

%results separating data by soz subjects (TLE vs Frontal):
%effects of any flicker on all patients:
results(end+1,:)={'all_soz-TLE',ied_results.soz_type.pct_est(strcmp(ied_results.soz_type.labels,'TLE')),ied_results.soz_type.pct_ci(strcmp(ied_results.soz_type.labels,'TLE'),:),ied_results.soz_type.pvals(strcmp(ied_results.soz_type.labels,'TLE')),ied_results.soz_type.df(strcmp(ied_results.soz_type.labels,'TLE'))};
results(end+1,:)=repmat({''},1,5);

%interaction of flicker modality and region of interest, for MTL SOZ subjects:
results(end+1,:)={'visual-flicker-in-visual_soz-TLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual'))};
results(end+1,:)={'visual-flicker-in-auditory_soz-TLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual'))};
results(end+1,:)={'visual-flicker-in-MTL_soz-TLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual'))};
results(end+1,:)={'visual-flicker-in-PFC_soz-TLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual'))};
results(end+1,:)=repmat({''},1,5);
results(end+1,:)={'audio-flicker-in-visual_soz-TLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio'))};
results(end+1,:)={'audio-flicker-in-auditory_soz-TLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio'))};
results(end+1,:)={'audio-flicker-in-MTL_soz-TLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio'))};
results(end+1,:)={'audio-flicker-in-PFC_soz-TLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio'))};
results(end+1,:)=repmat({''},1,5);
results(end+1,:)={'audiovisual-flicker-in-visual_soz-TLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV'))};
results(end+1,:)={'audiovisual-flicker-in-auditory_soz-TLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV'))};
results(end+1,:)={'audiovisual-flicker-in-MTL_soz-TLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV'))};
results(end+1,:)={'audiovisual-flicker-in-PFC_soz-TLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'TLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV'))};
results(end+1,:)=repmat({''},1,5);

%interaction of flicker modality and region of interest, for Frontal SOZ subjects:
results(end+1,:)={'all_soz-FLE',ied_results.soz_type.pct_est(strcmp(ied_results.soz_type.labels,'FLE')),ied_results.soz_type.pct_ci(strcmp(ied_results.soz_type.labels,'FLE'),:),ied_results.soz_type.pvals(strcmp(ied_results.soz_type.labels,'FLE')),ied_results.soz_type.df(strcmp(ied_results.soz_type.labels,'FLE'))};
results(end+1,:)=repmat({''},1,5);
results(end+1,:)={'visual-flicker-in-visual_soz-FLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual'))};
results(end+1,:)={'visual-flicker-in-auditory_soz-FLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual'))};
results(end+1,:)={'visual-flicker-in-MTL_soz-FLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual'))};
results(end+1,:)={'visual-flicker-in-PFC_soz-FLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Visual'))};
results(end+1,:)=repmat({''},1,5);
results(end+1,:)={'audio-flicker-in-visual_soz-FLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio'))};
results(end+1,:)={'audio-flicker-in-auditory_soz-FLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio'))};
results(end+1,:)={'audio-flicker-in-MTL_soz-FLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio'))};
results(end+1,:)={'audio-flicker-in-PFC_soz-FLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'Audio'))};
results(end+1,:)=repmat({''},1,5);
results(end+1,:)={'audiovisual-flicker-in-visual_soz-FLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'visual_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV'))};
results(end+1,:)={'audiovisual-flicker-in-auditory_soz-FLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'audio_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV'))};
results(end+1,:)={'audiovisual-flicker-in-MTL_soz-FLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'MTL_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV'))};
results(end+1,:)={'audiovisual-flicker-in-PFC_soz-FLE',ied_results.region_modality_soztle.pct_est(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')),[ied_results.region_modality_soztle.pct_ci_lower(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')) ied_results.region_modality_soztle.pct_ci_upper(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV'))],ied_results.region_modality_soztle.pvals(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV')),ied_results.region_modality_soztle.df(strcmp(ied_results.region_modality_soztle.soz_types,'FLE'),strcmp(ied_results.region_modality_soztle.regions,'PFC_regions'),strcmp(ied_results.region_modality_soztle.modalities,'AV'))};
results(end+1,:)=repmat({''},1,5);

%% find number of sessions and subjects included in TLE and FLE analyses (Figure_S9D):

metadata_tbl=readtable([root_dir '/FlickerStudyMetadata.xlsx'],'Sheet','Subjects','PreserveVariableNames',1);
[~,all_sessions]=fetch_flicker_subjectIDs(root_dir,'all');
all_sessions(strcmp(all_sessions.task,'spep'),:)=[];

TLE_sessions=all_sessions(ismember(all_sessions.sub,metadata_tbl.Subject_ID(strcmp(metadata_tbl.('Seizure focus broader classification'),'TLE'))),:);
FLE_sessions=all_sessions(ismember(all_sessions.sub,metadata_tbl.Subject_ID(strcmp(metadata_tbl.('Seizure focus broader classification'),'FLE'))),:);

disp(['Total number of TLE subjects: ' num2str(length(unique(TLE_sessions.sub)))]);
disp(['Total number of sessions from TLE subjects: ' num2str(size(TLE_sessions,1))]);

disp(['Total number of FLE subjects: ' num2str(length(unique(FLE_sessions.sub)))]);
disp(['Total number of sessions from FLE subjects: ' num2str(size(FLE_sessions,1))]);

%% save results table and ied table:

copyfile([root_dir '/stg-analyses/ied-analysis/IED_count_table_all-sessions.csv'],[source_data '/IED_count_table_all-sessions.csv']);
results_tbl=array2table([repmat({''},1,size(results,2));results],'VariableNames',{'analysis','ied_change','confidence_interval','p-value','df'});
temp_mat=[];
for i=1:length(results_tbl.confidence_interval)
    if isempty(results_tbl.confidence_interval{i})
        temp_mat=[temp_mat;nan nan];
    else
        temp_mat=[temp_mat;results_tbl.confidence_interval{i}];
    end
end
results_tbl=[results_tbl(:,1:2) array2table(num2cell(temp_mat),'VariableNames',{'confidence_interval_min','confidence_interval_max'}) results_tbl(:,4:5)];
results_tbl.analysis(arrayfun(@(x) isempty(x{:}),results_tbl.analysis))={'Figure_6B_1','Figure_6B_3','','','Figure_6B_2','Figure_S9B','Figure_S9C','','Figure_S9D_1','','','','Figure_S9D_2','','','',''}';
results_tbl(arrayfun(@(x) isempty(x{:}),results_tbl.analysis),:)=[];
writetable(results_tbl,[source_data '/Figure_6B_S9.csv']);

%% plot data of interest:

%flicker stim vs baseline, whole brain level:
figure_making('width',1,'height',2);
set(gcf,'color','w');
yline(0,'--');
hold on;
i=1;
if ~isempty(results{i,1})
    scatter(i,results{i,2},50,'k','filled');
    plot([i i],[results{i,3}(1) results{i,3}(2)],'k','LineWidth',1);
    
    asterisk=pval_to_asterisk(results{i,4});
    if ~isempty(asterisk)
        text(i-0.2,results{i,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
    end
end
ylabel('% IED rate change relative to baseline');
ax=gca;
ax.XAxis.Visible='off';
xlabel('Any flicker');
ax.XLabel.Visible='on';

figure_making('operation','save','filename',[figures_dir '/Figure_6B_1.pdf']);

%plot scatter of session averages:
%NEED ADDPATH FOR MATLAB VIOLIN PLOT FUNCTION REPO
figure_making('width',1.5,'height',2);
set(gcf,'color','w');
yline(0,'--');
sess_ied_values=ied_by_sess.all_chan.pct;
v=violinplot(sess_ied_values);
v.ViolinAlpha=0;
v.ViolinPlot.EdgeAlpha=0;
v.ScatterPlot.MarkerFaceAlpha=0.5;
v.ScatterPlot.MarkerFaceColor='flat';
v.ScatterPlot.CData=repmat([0 0 0],[length(sess_ied_values),1]);
v.ScatterPlot.SizeData=5;
v.BoxPlot.FaceAlpha=0;
v.BoxPlot.EdgeAlpha=0;
v.MedianPlot.Marker='none';
v.WhiskerPlot.Visible=0;
ylabel('% IED rate change relative to baseline');
ax=gca;
ax.XAxis.Visible='off';
xlabel('Any flicker');
ax.XLabel.Visible='on';
ylim([-25 25]);

%save figure and associated data:
figure_making('operation','save','filename',[figures_dir '/Figure_S9A.pdf']);
writetable(cell2table(num2cell(sess_ied_values),'VariableNames',{'percent_ied_rate_change_rel_baseline'}),[source_data '/Figure_S9A.csv']);

%flicker stim vs baseline, whole brain, modulated vs non-modulated:
figure_making('width',2,'height',2);
set(gcf,'color','w');
yline(0,'--');
hold on;
to_plot=[18,19];
for i=to_plot
    if ~isempty(results{i,1})
        scatter(i,results{i,2},50,'k','filled');
        plot([i i],[results{i,3}(1) results{i,3}(2)],'k','LineWidth',1);
        
        asterisk=pval_to_asterisk(results{i,4});
        if ~isempty(asterisk)
            text(i-0.2,results{i,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
        end
    end
end
ylabel('% IED rate change relative to baseline');
xlim([to_plot(1)-1 to_plot(2)+1]);
ax=gca;
ax.XTick=to_plot;
ax.XTickLabel={sprintf('non-mod/\\newlinelow-mod'),'high-mod'};
xlabel('Channels');

figure_making('operation','save','filename',[figures_dir '/Figure_6B_2.pdf']);


%interaction between modality and region:
figure_making('width',3,'height',2);
set(gcf,'color','w');
yline(0,'--');
hold on;
for i=3:6
    scatter(i-0.1,results{i,2},50,condition_color('V'),'filled');
    plot([i-0.1 i-0.1],[results{i,3}(1) results{i,3}(2)],'Color',condition_color('V'),'LineWidth',1);
    asterisk=pval_to_asterisk(results{i,4});
    if ~isempty(asterisk)
        text(i-0.3,results{i,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
    end
    
    scatter(i,results{i+10,2},50,condition_color('AV'),'filled');
    plot([i i],[results{i+10,3}(1) results{i+10,3}(2)],'Color',condition_color('AV'),'LineWidth',1);
    asterisk=pval_to_asterisk(results{i+10,4});
    if ~isempty(asterisk)
        text(i-0.2,results{i+10,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
    end
    
    scatter(i+0.1,results{i+5,2},50,condition_color('A'),'filled');
    plot([i+0.1 i+0.1],[results{i+5,3}(1) results{i+5,3}(2)],'Color',condition_color('A'),'LineWidth',1);
    asterisk=pval_to_asterisk(results{i+5,4});
    if ~isempty(asterisk)
        text(i-0.1,results{i+5,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
    end
end
ylabel('% IED rate change relative to baseline');
xlim([2 7]);
ylim([-60 60]);
ax=gca;
ax.XTick=3:6;
ax.XTickLabel={'Visual','Audio','MTL','PFC'};
xlabel('Region');
text(5.8,ax.YLim(2)-5,['\color[rgb]{' num2str(condition_color('V')) '}V  \color[rgb]{' num2str(condition_color('AV')) '}AV  \color[rgb]{' num2str(condition_color('A')) '}A'],'FontSize',7);

figure_making('operation','save','filename',[figures_dir '/Figure_6B_3.pdf']);


%interaction between modality and region for TLE SOZ subjects:
figure_making('width',3,'height',2);
set(gcf,'color','w');
yline(0,'--');
hold on;
i=90;
scatter(i+1,results{i,2},50,'k','filled');
plot([i+1 i+1],[results{i,3}(1) results{i,3}(2)],'k','LineWidth',1);
asterisk=pval_to_asterisk(results{i,4});
if ~isempty(asterisk)
    text(i+1-0.2,results{i,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
end
for i=92:95
    scatter(i-0.1,results{i,2},50,condition_color('V'),'filled');
    plot([i-0.1 i-0.1],[results{i,3}(1) results{i,3}(2)],'Color',condition_color('V'),'LineWidth',1);
    asterisk=pval_to_asterisk(results{i,4});
    if ~isempty(asterisk)
        text(i-0.3,results{i,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
    end
    
    scatter(i,results{i+10,2},50,condition_color('AV'),'filled');
    plot([i i],[results{i+10,3}(1) results{i+10,3}(2)],'Color',condition_color('AV'),'LineWidth',1);
    asterisk=pval_to_asterisk(results{i+10,4});
    if ~isempty(asterisk)
        if i==92
            text(i+0.35,results{i+10,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
        else
            text(i-0.2,results{i+10,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
        end
    end
    
    scatter(i+0.1,results{i+5,2},50,condition_color('A'),'filled');
    plot([i+0.1 i+0.1],[results{i+5,3}(1) results{i+5,3}(2)],'Color',condition_color('A'),'LineWidth',1);
    asterisk=pval_to_asterisk(results{i+5,4});
    if ~isempty(asterisk)
        if i==94 %so can see asterix better
            text(i+0.45,results{i+5,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
        else
            text(i-0.1,results{i+5,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
        end
    end
end
ylabel('% IED rate change relative to baseline');
xlim([90 96]);
ylim([-60 40]);
ax=gca;
ax.XTick=91:95;
ax.XTickLabel={'Whole brain','Visual','Audio','MTL','PFC'};
xlabel('Region');
text(93.8,ax.YLim(2),['\color[rgb]{' num2str([0 0 0]) '}any  \color[rgb]{' num2str(condition_color('V')) '}V  \color[rgb]{' num2str(condition_color('AV')) '}AV  \color[rgb]{' num2str(condition_color('A')) '}A'],'FontSize',7);

figure_making('operation','save','filename',[figures_dir '/Figure_S9D_1.pdf']);

%interaction between modality and region for FLE SOZ subjects:
figure_making('width',3,'height',2);
set(gcf,'color','w');
yline(0,'--');
hold on;
i=107;
scatter(i+1,results{i,2},50,'k','filled');
plot([i+1 i+1],[results{i,3}(1) results{i,3}(2)],'k','LineWidth',1);
asterisk=pval_to_asterisk(results{i,4});
if ~isempty(asterisk)
    text(i+1-0.2,results{i,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
end
for i=109:112
    scatter(i-0.1,results{i,2},50,condition_color('V'),'filled');
    plot([i-0.1 i-0.1],[results{i,3}(1) results{i,3}(2)],'Color',condition_color('V'),'LineWidth',1);
    asterisk=pval_to_asterisk(results{i,4});
    if ~isempty(asterisk)
        text(i-0.3,results{i,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
    end
    
    scatter(i,results{i+10,2},50,condition_color('AV'),'filled');
    plot([i i],[results{i+10,3}(1) results{i+10,3}(2)],'Color',condition_color('AV'),'LineWidth',1);
    asterisk=pval_to_asterisk(results{i+10,4});
    if ~isempty(asterisk)
        text(i-0.2,results{i+10,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
    end
    
    scatter(i+0.1,results{i+5,2},50,condition_color('A'),'filled');
    plot([i+0.1 i+0.1],[results{i+5,3}(1) results{i+5,3}(2)],'Color',condition_color('A'),'LineWidth',1);
    asterisk=pval_to_asterisk(results{i+5,4});
    if ~isempty(asterisk)
        if i==112
            text(i+0.45,results{i+5,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
        else
            text(i-0.1,results{i+5,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
        end
    end
end
ylabel('% IED rate change relative to baseline');
xlim([107 113]);
ylim([-50 100]);
ax=gca;
ax.XTick=108:112;
ax.XTickLabel={'Whole-brain','Visual','Audio','MTL','PFC'};
xlabel('Region');
text(110.8,ax.YLim(2),['\color[rgb]{' num2str([0 0 0]) '}any  \color[rgb]{' num2str(condition_color('V')) '}V  \color[rgb]{' num2str(condition_color('AV')) '}AV  \color[rgb]{' num2str(condition_color('A')) '}A'],'FontSize',7);

figure_making('operation','save','filename',[figures_dir '/Figure_S9D_2.pdf']);


%showing by flickerneuro condition:
figure_making('width',5.5,'height',2);
set(gcf,'color','w');
yline(0,'--');
hold on;
for i=21:32
    if ~isempty(results{i,1})
        temp=strsplit(results{i,1},'_');
        scatter(i,results{i,2},50,condition_color(temp{2}),'filled');
        plot([i i],[results{i,3}(1) results{i,3}(2)],'Color',condition_color(temp{2}),'LineWidth',1);
        asterisk=pval_to_asterisk(results{i,4});
        if ~isempty(asterisk)
            text(i-0.2,results{i,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
        end
    end
    if mod(i-21,3)==2 && i~=32
        xline(i+0.5,'-');
    end
end
ylabel('% IED rate change relative to baseline');
xlim([20 33]);
ax=gca;
ax.XTick=[22,25,28,31];
ax.XTickLabel={'5.5Hz','40Hz','80Hz','Random'};
xlabel('Stimulation frequency');
text(29.9,ax.YLim(2)-1,['\color[rgb]{' num2str(condition_color('V')) '}V         \color[rgb]{' num2str(condition_color('AV')) '}AV        \color[rgb]{' num2str(condition_color('A')) '}A'],'FontSize',7);

figure_making('operation','save','filename',[figures_dir '/Figure_S9B.pdf']);


%flickerfreq results:
figure_making('width',7.5,'height',2);
set(gcf,'color','w');
yline(0,'--');
hold on;
for i=34:60
    if ~isempty(results{i,1})
        temp=strsplit(results{i,1},'_');
        scatter(i-0.1,results{i,2},50,condition_color(temp{2}),'filled');
        plot([i-0.1 i-0.1],[results{i,3}(1) results{i,3}(2)],'Color',condition_color(temp{2}),'LineWidth',1);
        asterisk=pval_to_asterisk(results{i,4});
        if ~isempty(asterisk)
            text(i-0.3,results{i,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
        end
    end
    
    if ~isempty(results{i,1})
        temp=strsplit(results{i+28,1},'_');
        scatter(i+0.1,results{i+28,2},50,condition_color(temp{2}),'filled');
        plot([i+0.1 i+0.1],[results{i+28,3}(1) results{i+28,3}(2)],'Color',condition_color(temp{2}),'LineWidth',1);
        asterisk=pval_to_asterisk(results{i+28,4});
        if ~isempty(asterisk)
            text(i-0.1,results{i+28,2},asterisk,'Rotation',90,'HorizontalAlignment','center','FontSize',12);
        end
    end
end
ylabel('% IED rate change relative to baseline');
xlim([33 61]);
ax=gca;
ax.XTick=34:60;
ax.XTickLabel=regexprep(flickerfreq_conditions,'Hz','');
xlabel('Stimulation frequency (Hz)');
text(59,ax.YLim(2)-5,['\color[rgb]{' num2str(condition_color('V')) '}V   \color[rgb]{' num2str(condition_color('A')) '}A'],'FontSize',7);

figure_making('operation','save','filename',[figures_dir '/Figure_S9C.pdf']);

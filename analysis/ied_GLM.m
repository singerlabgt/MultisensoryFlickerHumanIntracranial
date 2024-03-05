


%% input and output directories:

opts=detectImportOptions([output_dir '/IED_count_table_all-sessions.csv']);
opts.VariableTypes(1)={'char'};
data_table_original=readtable([output_dir '/IED_count_table_all-sessions.csv'],opts);


%% Generate a data table where ied counts are aggregated
% for channels from same consecutive condition;
% used for analysis examining variability across individual channels
% (e.g. IED change vs. flicker modulation strength)


gen_chan_table = true;

data_table=data_table_original;

subject_chans = unique(data_table(:,[1 5]),'rows');

ied_table_chan = [];  %table with IED counts per channel, summed by trial

disp('Generating channel table (may take a few minutes)...')

for kk = 1:size(subject_chans,1)
    sel_inds = and(strcmp(data_table.subjectID,subject_chans.subjectID{kk}),strcmp(data_table.ch,subject_chans.ch{kk}));
    chan_table = data_table(sel_inds,:);
    
    conds_temp = unique(chan_table.cond);
    for kl = 1:length(conds_temp)
        chan_table_cond = chan_table(strcmp(chan_table.cond,conds_temp{kl}),:);
        baseline_inds = logical(chan_table_cond.baseline);
        
        ied_count_baseline = sum(chan_table_cond.ied_count(baseline_inds));
        ied_count_flicker = sum(chan_table_cond.ied_count(~baseline_inds));
        
        temprow = chan_table_cond(1,:);
        temprow.trial_nber = 'all'; temprow.baseline = 1; temprow.ied_count = ied_count_baseline;
        ied_table_chan = [ied_table_chan; temprow];
        
        temprow.trial_nber = 'all'; temprow.baseline = 0; temprow.ied_count = ied_count_flicker;
        ied_table_chan = [ied_table_chan; temprow];
        
        
    end
end
ied_change = zeros(size(ied_table_chan.ied_count)); 
ied_change_temp = ied_table_chan.ied_count(~logical(ied_table_chan.baseline)) - ied_table_chan.ied_count(logical(ied_table_chan.baseline));
ied_change(1:2:length(ied_change)) = ied_change_temp;
ied_change(2:2:length(ied_change)) = ied_change_temp;
ied_table_chan.ied_change = ied_change;

save([output_dir '/ied_table_chan.mat'],'ied_table_chan');

%% Create whole brain data table, where IED counts are summed
% across channels for different sets of brain regions.
    %prefrontal, mesial temporal, audio, visual, other
% 
% used for most IED analysis, e.g. examining whether IED reduction was
% significant for different regions and/or modalities

disp('Generating whole-brain IED tables...')

for file_type_index = 1:3
    data_table=data_table_original;
    if file_type_index == 1
        gen_wb_table = true;
        task2only = false; 
        task1only = false;
    elseif file_type_index == 2
        task1only = true;
        task2only = false; 
    else
        task2only = true; 
        task1only = false; 
    end

%task 2 = sweep of all frequencies at shorter duration; task 1 = original
%smaller frequency grid
if task1only
    data_table(strcmp(data_table.task_name,'flickerfreq'),:) = [];
elseif task2only
    data_table(strcmp(data_table.task_name,'flickerneuro'),:) = [];
else
    
end

if gen_wb_table
    %for old table (v2) use: subject_sessions = unique(data_table(:,[1 2]),'rows');
    subject_sessions = unique(data_table(:,[1 3]),'rows');
    
    ied_table_wb = table();   %table with IED counts summed across all chans
    
    n_rows = 10856;
    
    ied_count_soz = zeros(n_rows,1);
    ied_count_nonsoz = zeros(n_rows,1);
    
    ied_count_pfc = zeros(n_rows,1);
    ied_count_mtl = zeros(n_rows,1);
    ied_count_audio = zeros(n_rows,1);
    ied_count_visual = zeros(n_rows,1);
    ied_count_other = zeros(n_rows,1);
    
    ied_counter = 1;
    for kk = 1:size(subject_sessions,1)
        sel_inds = and(strcmp(data_table.subjectID,subject_sessions.subjectID{kk}),data_table.ses_nber == subject_sessions.ses_nber(kk));
        
        sel_chans = unique(data_table.ch(sel_inds));
        
        ied_baseline_temp = 0;
        ied_flicker_temp = 0;
        
        sub_table = data_table(sel_inds,:);
        cond_temp = unique(sub_table.cond);
        
        if task1only && length(cond_temp)>15  %skip task 2 with more flicker frequencies
            continue
        end
        
        for kl = 1:length(sel_chans)
            sub_table_chan = sub_table(strcmp(sub_table.ch,sel_chans{kl}),:);
            
            if kl == 1
                ied_inds = ied_counter:ied_counter + size(sub_table_chan,1) -1;
                ied_count_temp = sub_table_chan.ied_count;
            else
                if size(sub_table_chan,1) == length(ied_count_temp)
                    ied_count_temp = ied_count_temp + sub_table_chan.ied_count;
                else
                    
                end
            end
            
            if any(strcmp(sub_table_chan.ch_region,'PFC_regions')) && size(sub_table_chan,1) == length(ied_count_temp)
                ied_count_pfc(ied_inds) = ied_count_pfc(ied_inds) + sub_table_chan.ied_count;
            end
            if any(strcmp(sub_table_chan.ch_region,'MTL_regions')) && size(sub_table_chan,1) == length(ied_count_temp)
                ied_count_mtl(ied_inds) = ied_count_mtl(ied_inds) + sub_table_chan.ied_count;
            end
            if any(strcmp(sub_table_chan.ch_region,'audio_regions')) && size(sub_table_chan,1) == length(ied_count_temp)
                ied_count_audio(ied_inds) = ied_count_audio(ied_inds) + sub_table_chan.ied_count;
            end
            if any(strcmp(sub_table_chan.ch_region,'visual_regions')) && size(sub_table_chan,1) == length(ied_count_temp)
                ied_count_visual(ied_inds) = ied_count_visual(ied_inds) + sub_table_chan.ied_count;
            end
            if any(strcmp(sub_table_chan.ch_region,'other')) && size(sub_table_chan,1) == length(ied_count_temp)
                ied_count_other(ied_inds) = ied_count_other(ied_inds) + sub_table_chan.ied_count;
            end
        end
        
        sub_table_chan.ied_count = ied_count_temp;
        for kl = 1:size(sub_table_chan,1)
            sub_table_chan.ch{kl} = 'all';
            sub_table_chan.ch_anat{kl} = 'all';
            sub_table_chan.ch_region{kl} = 'all';
        end
        ied_table_wb = [ied_table_wb; sub_table_chan];
        ied_counter = ied_counter + size(sub_table_chan,1);
    end
    
    
    %ied_table_wb.ied_change = ied_change;
    if task1only || task2only
        n_rows = size(ied_table_wb,1);
        ied_count_soz = ied_count_soz(1:n_rows);
        ied_count_pfc = ied_count_pfc(1:n_rows);
        ied_count_mtl = ied_count_mtl(1:n_rows);
        ied_count_visual = ied_count_visual(1:n_rows);
        ied_count_audio = ied_count_audio(1:n_rows);
        ied_count_other = ied_count_other(1:n_rows);
    end
    
    ied_table_wb.ied_count_soz = ied_count_soz;
    ied_table_wb.ied_count_nonsoz = ied_table_wb.ied_count-ied_count_soz;
    ied_table_wb.ied_count_pfc = ied_count_pfc;
    ied_table_wb.ied_count_mtl = ied_count_mtl;
    ied_table_wb.ied_count_visual = ied_count_visual;
    ied_table_wb.ied_count_audio = ied_count_audio;
    ied_table_wb.ied_count_other = ied_count_other;
    
    ied_change = zeros(n_rows,1);
    ied_change_temp = ied_table_wb.ied_count(~logical(ied_table_wb.baseline)) - ied_table_wb.ied_count(logical(ied_table_wb.baseline));
    ied_change(1:2:length(ied_change)) = ied_change_temp;
    ied_change(2:2:length(ied_change)) = ied_change_temp;
    
    ied_table_wb.modality = get_modalities(ied_table_wb.cond);
    
    if task1only
        save([output_dir, '/IED_table_wb_task1.mat'],'ied_table_wb')
    elseif task2only
        save([output_dir, '/IED_table_wb_task2.mat'],'ied_table_wb')
    else
        save([output_dir, '/IED_table_wb.mat'],'ied_table_wb')
    end

end

end

%% initialize structure to append results for different analysis (will be saved later)

disp('Starting stats:')

ied_results = struct();

%% Stats: whole-brain IED count vs. any flicker (Figure_6b_1, all patients + by SOZ)

load([output_dir, '/IED_table_wb.mat'])

soz_types = unique(ied_table_wb.soz_anat);

%soz_types = {'All'};
soz_types = [{'All'}; soz_types];

soz_est = zeros(size(soz_types));
soz_se = zeros(size(soz_types));
soz_pvals = zeros(size(soz_types));
soz_baseline = zeros(size(soz_types));
soz_df = zeros(size(soz_types));
soz_ci = zeros(length(soz_types),2);

for kk = 1:length(soz_types)
    if strcmp(soz_types{kk},'All')
        sub_table = ied_table_wb;
    else
        sub_table = ied_table_wb;
        pts_to_exc = unique(sub_table.subjectID(~strcmp(sub_table.soz_anat,soz_types{kk})));
        sub_table(contains(sub_table.subjectID,pts_to_exc),:) = [];
    end
    %sub_table.baseline = 1 - sub_table.baseline;
    glme = fitglme(sub_table,'ied_count ~ 1 + baseline + (1|subjectID) + (1|ses_nber)', ...
        'Distribution','Poisson','Link','log');
    display(glme)
    
    soz_pvals(kk) = glme.Coefficients.pValue(2);
    soz_est(kk) = glme.Coefficients.Estimate(2);
    soz_se(kk) = glme.Coefficients.SE(2);
    soz_baseline(kk) = exp(glme.Coefficients.Estimate(1));
    soz_df(kk) = glme.Coefficients.DF(2);
    
    soz_ci(kk,1) = glme.Coefficients.Lower(2);
    soz_ci(kk,2) = glme.Coefficients.Upper(2);
end

ied_mean = -100*(exp(soz_est)-1);
ied_serr = -100*(exp(soz_se)-1);


temp_struct = struct();
temp_struct.labels = soz_types;

temp_struct.pvals = soz_pvals;
temp_struct.est = soz_est;
temp_struct.se = soz_se;
temp_struct.ci = soz_ci;
temp_struct.baseline = soz_baseline;
temp_struct.df = soz_df;

temp_struct.pct_est= -100*(exp(soz_est)-1);
temp_struct.pct_se = -100*(exp(soz_se)-1);
temp_struct.pct_ci = -100*(exp(soz_ci)-1);

ied_results.soz_type = temp_struct;

%% Stats: does stronger flicker modulation predict greater IED change? Figure_6B_2

load([output_dir, '/IED_table_chan.mat'])

mod_types = {'Low mod','High mod'};

ied_table_chan.mod_high = zeros(size(ied_table_chan,1),1);

mod_thresh_high = 2;
mod_thresh_low = -10000000;

unique_chans = unique(ied_table_chan(:,[1,5]),'rows');
for kk = size(unique_chans,1):-1:1
    sub_inds = and(strcmp(ied_table_chan.subjectID,unique_chans.subjectID{kk}),strcmp(ied_table_chan.ch,unique_chans.ch{kk}));
    sub_table = ied_table_chan(sub_inds,:);
    
    avg_mod = nanmean(sub_table.mod_amp);
    if avg_mod< mod_thresh_low
        ied_table_chan(sub_inds,:) = [];
        
    elseif avg_mod< mod_thresh_high
        ied_table_chan.mod_amp(sub_inds) = avg_mod;
        
    else
        ied_table_chan.mod_amp(sub_inds) = avg_mod;
        ied_table_chan.mod_high(sub_inds) = 1;
    end
end

pvals = zeros(size(mod_types));
ests = zeros(size(mod_types));
est_se = zeros(size(mod_types));
est_ci = zeros(length(mod_types),2);
baseline = zeros(size(mod_types));
df = zeros(size(mod_types));

n_chan_low = length(unique(ied_table_chan.ch(ied_table_chan.mod_high == 0)));
n_chan_high = length(unique(ied_table_chan.ch(ied_table_chan.mod_high == 1)));

new_thresh = 1.5;
ied_table_chan.mod_high(ied_table_chan.mod_amp<new_thresh,:) = 0;
ied_table_chan.mod_high(ied_table_chan.mod_amp>new_thresh,:) = 1;

%ied_table_chan(ied_table_chan.mod_sig>0.05,:) = [];

glme = fitglme(ied_table_chan(ied_table_chan.mod_high==0,:), 'ied_count ~ 1 + baseline + (1|subjectID)',...
                    'Distribution','Poisson','Link','log');
pvals(1) = glme.Coefficients(2,6).pValue;
ests(1) = glme.Coefficients.Estimate(2);
est_se(1) =  glme.Coefficients.SE(2);
baseline(1) = exp(glme.Coefficients.Estimate(1));
est_ci(1,:) = [glme.Coefficients.Lower(2),glme.Coefficients.Upper(2)];
df(1) = glme.Coefficients.DF(2);
display(glme)

glme = fitglme(ied_table_chan(ied_table_chan.mod_high==1,:), 'ied_count ~ 1 + baseline + (1|subjectID)',...
                    'Distribution','Poisson','Link','log');
display(glme)

ests(2) = glme.Coefficients.Estimate(2);
est_se(2) =  glme.Coefficients.SE(2);
baseline(2) = exp(glme.Coefficients.Estimate(1));
pvals(2) = glme.Coefficients(2,6).pValue;
est_ci(2,:) = [glme.Coefficients.Lower(2),glme.Coefficients.Upper(2)];
df(2) = glme.Coefficients.DF(2);

glme = fitglme(ied_table_chan, 'ied_count ~ 1 + baseline*mod_high + (1|subjectID)',...
                    'Distribution','Poisson','Link','log');
                
display(glme)

p_diff = glme.Coefficients.pValue(4);
df_diff =  glme.Coefficients.DF(4);

temp_struct = struct();
temp_struct.labels = mod_types;
 
temp_struct.pvals = pvals;
temp_struct.p_compare = p_diff;
temp_struct.df_compare = df_diff;
temp_struct.est = ests;
temp_struct.se = est_se;
temp_struct.ci = est_ci;
temp_struct.baseline = exp(baseline);
temp_struct.df = df;
temp_struct.df_diff = df_diff;

temp_struct.pct_est= -100*(exp(ests)-1);
temp_struct.pct_se = -100*(exp(est_se)-1);
temp_struct.pct_ci = -100*(exp(est_ci)-1);

ied_results.flicker_strength = temp_struct;


%% Interaction of modality and IED change in diff. regions. Figure_6B_3
load([output_dir, '/IED_table_wb.mat'])
soz_types = {'MTL_regions','PFC_regions','visual_regions','audio_regions','other'};
soz_types_plot = {'MTL','PFC','Visual','Audio','Other'};
mod_types = {'Visual','Audio','AV'};

ied_count_baseline = cell(length(soz_types),length(mod_types));
ied_count_flicker = cell(length(soz_types),length(mod_types));
ied_pct_change = cell(length(soz_types),length(mod_types));

ied_pct_by_sess = cell(length(soz_types),length(mod_types));
ied_base_by_sess = cell(length(soz_types),length(mod_types));
ied_fl_by_sess = cell(length(soz_types),length(mod_types));

modal_est = zeros(length(soz_types),length(mod_types));
modal_se = zeros(length(soz_types),length(mod_types));
modal_baseline = zeros(length(soz_types),length(mod_types));
modal_pvals = zeros(length(soz_types),length(mod_types));
modal_ci_upper = zeros(length(soz_types),length(mod_types));
modal_ci_lower = zeros(length(soz_types),length(mod_types));
modal_df = zeros(length(soz_types),length(mod_types));

subs = unique(ied_table_wb.subjectID);

ied_means_pfc = nan(3,2);
ied_means_visual = nan(3,2);
ied_means_audio = nan(3,2);

for kk = 1:length(soz_types)
    for kk2 = 1:length(mod_types)
        ied_table_temp = ied_table_wb(strcmp(ied_table_wb.modality,mod_types{kk2}),:);
        if kk == 1
            glme = fitglme(ied_table_temp, 'ied_count_mtl ~ 1 + baseline + (1|subjectID)',...
                    'Distribution','Poisson','Link','log');
            modal_pvals(kk,kk2) = glme.Coefficients.pValue(2);
            modal_est(kk,kk2) = glme.Coefficients.Estimate(2);
            modal_se(kk,kk2) = glme.Coefficients.SE(2);
            modal_ci_lower(kk,kk2) = glme.Coefficients.Lower(2);
            modal_ci_upper(kk,kk2) = glme.Coefficients.Upper(2);
            modal_baseline(kk,kk2) = exp(glme.Coefficients.Estimate(1));
            modal_df(kk,kk2) = glme.Coefficients.DF(2);
            
        elseif kk == 2
            glme = fitglme(ied_table_temp, 'ied_count_pfc ~ 1 + baseline + (1|subjectID)',...
                    'Distribution','Poisson','Link','log');
            modal_pvals(kk,kk2) = glme.Coefficients.pValue(2);
            modal_est(kk,kk2) = glme.Coefficients.Estimate(2);
            modal_se(kk,kk2) = glme.Coefficients.SE(2);
            modal_ci_lower(kk,kk2) = glme.Coefficients.Lower(2);
            modal_ci_upper(kk,kk2) = glme.Coefficients.Upper(2);
            modal_baseline(kk,kk2) = exp(glme.Coefficients.Estimate(1));
            modal_df(kk,kk2) = glme.Coefficients.DF(2);
            
            ied_means_pfc(kk2,1) = mean(ied_table_temp.ied_count_pfc(ied_table_temp.baseline == 0));
            ied_means_pfc(kk2,2) = mean(ied_table_temp.ied_count_pfc(ied_table_temp.baseline == 1));
        elseif kk==4
            glme = fitglme(ied_table_temp, 'ied_count_audio ~ 1 + baseline + (1|subjectID)',...
                    'Distribution','Poisson','Link','log');
            modal_pvals(kk,kk2) = glme.Coefficients.pValue(2);
            modal_est(kk,kk2) = glme.Coefficients.Estimate(2);
            modal_se(kk,kk2) = glme.Coefficients.SE(2);
            modal_ci_lower(kk,kk2) = glme.Coefficients.Lower(2);
            modal_ci_upper(kk,kk2) = glme.Coefficients.Upper(2);
            modal_baseline(kk,kk2) = exp(glme.Coefficients.Estimate(1));
            modal_df(kk,kk2) = glme.Coefficients.DF(2);
            
            ied_means_audio(kk2,1) = mean(ied_table_temp.ied_count_audio(ied_table_temp.baseline == 0));
            ied_means_audio(kk2,2) = mean(ied_table_temp.ied_count_audio(ied_table_temp.baseline == 1));
        elseif kk==3
            glme = fitglme(ied_table_temp, 'ied_count_visual ~ 1 + baseline + (1|subjectID)',...
                    'Distribution','Poisson','Link','log');
            modal_pvals(kk,kk2) = glme.Coefficients.pValue(2);
            modal_est(kk,kk2) = glme.Coefficients.Estimate(2);
            modal_se(kk,kk2) = glme.Coefficients.SE(2);
            modal_ci_lower(kk,kk2) = glme.Coefficients.Lower(2);
            modal_ci_upper(kk,kk2) = glme.Coefficients.Upper(2);
            modal_baseline(kk,kk2) = exp(glme.Coefficients.Estimate(1));
            modal_df(kk,kk2) = glme.Coefficients.DF(2);
            
            ied_means_visual(kk2,1) = mean(ied_table_temp.ied_count_visual(ied_table_temp.baseline == 0));
            ied_means_visual(kk2,2) = mean(ied_table_temp.ied_count_visual(ied_table_temp.baseline == 1));
        elseif kk==5
            glme = fitglme(ied_table_temp, 'ied_count_other ~ 1 + baseline + (1|subjectID)',...
                    'Distribution','Poisson','Link','log');
            modal_pvals(kk,kk2) = glme.Coefficients.pValue(2);
            modal_est(kk,kk2) = glme.Coefficients.Estimate(2);
            modal_se(kk,kk2) = glme.Coefficients.SE(2);
            modal_ci_lower(kk,kk2) = glme.Coefficients.Lower(2);
            modal_ci_upper(kk,kk2) = glme.Coefficients.Upper(2);
            modal_baseline(kk,kk2) = exp(glme.Coefficients.Estimate(1));
            modal_df(kk,kk2) = glme.Coefficients.DF(2);
        end
        for kl = 1:length(subs)
            sub_table = ied_table_wb(strcmp(ied_table_wb.subjectID,subs{kl}),:);
            sub_table = sub_table(strcmp(sub_table.modality,mod_types{kk2}),:);
            if kk == 1
                ied_count_baseline{kk,kk2} = [ied_count_baseline{kk,kk2}; mean(sub_table.ied_count_mtl(sub_table.baseline == 1))];
                ied_count_flicker{kk,kk2} = [ied_count_flicker{kk,kk2}; mean(sub_table.ied_count_mtl(sub_table.baseline == 0))];
            elseif kk == 2
                ied_count_baseline{kk,kk2} = [ied_count_baseline{kk,kk2}; mean(sub_table.ied_count_pfc(sub_table.baseline == 1))];
                ied_count_flicker{kk,kk2} = [ied_count_flicker{kk,kk2}; mean(sub_table.ied_count_pfc(sub_table.baseline == 0))];
            elseif kk==3
                ied_count_baseline{kk,kk2} = [ied_count_baseline{kk,kk2}; mean(sub_table.ied_count_audio(sub_table.baseline == 1))];
                ied_count_flicker{kk,kk2} = [ied_count_flicker{kk,kk2}; mean(sub_table.ied_count_audio(sub_table.baseline == 0))];
            elseif kk==4
                ied_count_baseline{kk,kk2} = [ied_count_baseline{kk,kk2}; mean(sub_table.ied_count_visual(sub_table.baseline == 1))];
                ied_count_flicker{kk,kk2} = [ied_count_flicker{kk,kk2}; mean(sub_table.ied_count_visual(sub_table.baseline == 0))];
            elseif kk==5
                ied_count_baseline{kk,kk2} = [ied_count_baseline{kk,kk2}; mean(sub_table.ied_count_other(sub_table.baseline == 1))];
                ied_count_flicker{kk,kk2} = [ied_count_flicker{kk,kk2}; mean(sub_table.ied_count_other(sub_table.baseline == 0))];
            end
        end
        
        unique_sess = unique(ied_table_temp(:,1:3),'rows');
        for kl = 1:size(unique_sess,1)
            tempinds = and(strcmp(ied_table_temp.subjectID,unique_sess.subjectID{kl}),strcmp(ied_table_temp.task_name,unique_sess.task_name{kl}));
            tempinds = and(tempinds,ied_table_temp.ses_nber == unique_sess.ses_nber(kk));
            subtable = ied_table_temp(tempinds,:);
            if kk == 1
                ied_pct_by_sess{kk,kk2} = [ied_pct_by_sess{kk,kk2}; 100*((nanmean(subtable.ied_count_mtl(subtable.baseline == 0))-nanmean(subtable.ied_count_mtl(subtable.baseline == 1)))./nanmean(subtable.ied_count_mtl(subtable.baseline == 1)))];
                ied_base_by_sess{kk,kk2} = [ied_base_by_sess{kk,kk2}; nanmean(subtable.ied_count_mtl(subtable.baseline == 1))];
                ied_fl_by_sess{kk,kk2} = [ied_fl_by_sess{kk,kk2};nanmean(subtable.ied_count_mtl(subtable.baseline == 0))];
            elseif kk == 2
                ied_pct_by_sess{kk,kk2} = [ied_pct_by_sess{kk,kk2}; 100*((nanmean(subtable.ied_count_pfc(subtable.baseline == 0))-nanmean(subtable.ied_count_pfc(subtable.baseline == 1)))./nanmean(subtable.ied_count_pfc(subtable.baseline == 1)))];
                ied_base_by_sess{kk,kk2} = [ied_base_by_sess{kk,kk2}; nanmean(subtable.ied_count_pfc(subtable.baseline == 1))];
                ied_fl_by_sess{kk,kk2} = [ied_fl_by_sess{kk,kk2}; nanmean(subtable.ied_count_pfc(subtable.baseline == 0))];
            elseif kk == 3
                ied_pct_by_sess{kk,kk2} = [ied_pct_by_sess{kk,kk2}; 100*((nanmean(subtable.ied_count_audio(subtable.baseline == 0))-nanmean(subtable.ied_count_audio(subtable.baseline == 1)))./nanmean(subtable.ied_count_audio(subtable.baseline == 1)))];
                ied_base_by_sess{kk,kk2} = [ied_base_by_sess{kk,kk2}; nanmean(subtable.ied_count_audio(subtable.baseline == 1))];
                ied_fl_by_sess{kk,kk2} = [ied_fl_by_sess{kk,kk2}; nanmean(subtable.ied_count_audio(subtable.baseline == 0))];
            elseif kk==4
                ied_pct_by_sess{kk,kk2} = [ied_pct_by_sess{kk,kk2}; 100*((nanmean(subtable.ied_count_visual(subtable.baseline == 0))-nanmean(subtable.ied_count_visual(subtable.baseline == 1)))./nanmean(subtable.ied_count_visual(subtable.baseline == 1)))];
                ied_base_by_sess{kk,kk2} = [ied_base_by_sess{kk,kk2}; nanmean(subtable.ied_count_visual(subtable.baseline == 1))];
                ied_fl_by_sess{kk,kk2} = [ied_fl_by_sess{kk,kk2}; nanmean(subtable.ied_count_visual(subtable.baseline == 0))];
            elseif kk == 5
                ied_pct_by_sess{kk,kk2} = [ied_pct_by_sess{kk,kk2}; 100*((nanmean(subtable.ied_count_other(subtable.baseline == 0))-nanmean(subtable.ied_count_other(subtable.baseline == 1)))./nanmean(subtable.ied_count_other(subtable.baseline == 1)))];
                ied_base_by_sess{kk,kk2} = [ied_base_by_sess{kk,kk2}; nanmean(subtable.ied_count_pfc(subtable.baseline == 1))];
                ied_fl_by_sess{kk,kk2} = [ied_fl_by_sess{kk,kk2}; nanmean(subtable.ied_count_pfc(subtable.baseline == 0))];
            end
        end
    end
end
for kk = 1:length(soz_types)
    for kl = 1:length(mod_types)
        ied_pct_change{kk,kl} = (ied_count_flicker{kk,kl} - ied_count_baseline{kk,kl});
    end
end

temp_struct = struct();
temp_struct.labels = soz_types_plot;
temp_struct.modalities = mod_types;
temp_struct.pvals = modal_pvals;
temp_struct.est = modal_est;
temp_struct.se = modal_se;
temp_struct.ci_upper = modal_ci_upper;
temp_struct.ci_lower = modal_ci_lower;
temp_struct.baseline = modal_baseline;
temp_struct.df = modal_df;

temp_struct.pct_est= -100*(exp(modal_est)-1);
temp_struct.pct_se = -100*(exp(modal_se)-1);
temp_struct.pct_ci_upper = -100*(exp(modal_ci_upper)-1);
temp_struct.pct_ci_lower = -100*(exp(modal_ci_lower)-1);

ied_results.region_modality = temp_struct;


ied_by_sess.region_modality.pct = ied_pct_by_sess;
ied_by_sess.region_modality.baseline = ied_base_by_sess;
ied_by_sess.region_modality.flicker = ied_fl_by_sess;

%% Interaction between regions and different flicker conditions (Supplement, B)
load([output_dir, '/IED_table_wb_task1.mat'])

unique_sessions = unique(ied_table_wb(:,[1 2 3]),'rows');
exc_inds = zeros(size(ied_table_wb,1),1);
for kk = size(unique_sessions,1):-1:1
    sub_inds = and(strcmp(ied_table_wb.subjectID,unique_sessions.subjectID{kk}),ied_table_wb.ses_nber==unique_sessions.ses_nber(kk));
    cond_temp = unique(ied_table_wb.cond(sub_inds));
    if length(cond_temp)>15
        exc_inds(sub_inds) = 1;
    end
end
ied_table_wb(logical(exc_inds),:) = [];

conditions = unique(ied_table_wb.cond);

pvals = zeros(size(conditions));
c_est = zeros(size(conditions));
c_se = zeros(size(conditions));
c_ci = zeros(length(conditions),2);
c_baseline = zeros(size(conditions));
c_df = zeros(size(conditions));

pct_est = zeros(size(conditions));
pct_se = zeros(size(conditions));
pct_ci = zeros(length(conditions),2);

for kk = 1:length(conditions)
    sub_table = ied_table_wb(strcmp(ied_table_wb.cond,conditions{kk}),:);
    
    glme = fitglme(sub_table, 'ied_count ~ 1 + baseline + (1|subjectID)',...
                    'Distribution','Poisson','Link','log');
                
    pvals(kk) = glme.Coefficients.pValue(2);
    c_est(kk) = glme.Coefficients.Estimate(2);
    c_baseline(kk) = exp(glme.Coefficients.Estimate(1));
    c_se(kk) =  glme.Coefficients.SE(2);
    c_ci(kk,1) = glme.Coefficients.Lower(2);
    c_ci(kk,2) = glme.Coefficients.Upper(2);
    c_df(kk) = glme.Coefficients.DF(2);
    
    pct_est(kk) = -100*(exp(glme.Coefficients.Estimate(2))-1);
    pct_se(kk) =  -100*(exp(glme.Coefficients.SE(2))-1);
    pct_ci(kk,1) = -100*(exp(glme.Coefficients.Lower(2))-1);
    pct_ci(kk,2) = -100*(exp(glme.Coefficients.Upper(2))-1);
end

temp_struct = struct();
temp_struct.labels = conditions;
temp_struct.pvals = pvals;
temp_struct.est = c_est;
temp_struct.se = c_se;
temp_struct.ci = c_ci;
temp_struct.baseline = c_baseline;
temp_struct.pct_est = pct_est;
temp_struct.pct_se = pct_se;
temp_struct.pct_ci = pct_ci;
temp_struct.df = c_df;

ied_results.condition_exp1 = temp_struct;

%% Same thing but now experiment 2 (full range of frequencies) (Supplement, C)
load([output_dir, '/IED_table_wb_task2.mat'])

conditions = unique(ied_table_wb.cond);

pvals = zeros(size(conditions));
c_est = zeros(size(conditions));
c_se = zeros(size(conditions));
c_ci = zeros(length(conditions),2);
c_baseline = zeros(size(conditions));
c_df = zeros(size(conditions));

pct_est = zeros(size(conditions));
pct_se = zeros(size(conditions));
pct_ci = zeros(length(conditions),2);

for kk = 1:length(conditions)
    sub_table = ied_table_wb(strcmp(ied_table_wb.cond,conditions{kk}),:);
    
    glme = fitglme(sub_table, 'ied_count ~ 1 + baseline + (1|subjectID)',...
                    'Distribution','Poisson','Link','log');
                
    pvals(kk) = glme.Coefficients.pValue(2);
    c_est(kk) = glme.Coefficients.Estimate(2);
    c_se(kk) =  glme.Coefficients.SE(2);
    c_df(kk) = glme.Coefficients.DF(2);
    c_ci(kk,1) = glme.Coefficients.Lower(2);
    c_ci(kk,2) = glme.Coefficients.Upper(2);
    
    pct_est(kk) = -100*(exp(glme.Coefficients.Estimate(2))-1);
    pct_se(kk) =  -100*(exp(glme.Coefficients.SE(2))-1);
    pct_ci(kk,1) = -100*(exp(glme.Coefficients.Lower(2))-1);
    pct_ci(kk,2) = -100*(exp(glme.Coefficients.Upper(2))-1);
end

temp_struct = struct();
temp_struct.labels = conditions;
temp_struct.pvals = pvals;
temp_struct.est = c_est;
temp_struct.se = c_se;
temp_struct.ci = c_ci;
temp_struct.df = c_df;
temp_struct.baseline = c_baseline;
temp_struct.pct_est = pct_est;
temp_struct.pct_se = pct_se;
temp_struct.pct_ci = pct_ci;

ied_results.condition_exp2 = temp_struct;

%% Interaction between region, flicker condition, and soz type (IED supplement, D)

load([output_dir, '/IED_table_wb.mat'])
soz_types = unique(ied_table_wb.soz_anat);

regions = {'MTL_regions','PFC_regions','visual_regions','audio_regions','other'};
regions_plot = {'MTL','PFC','Visual','Audio','Other'};
mod_types = {'Visual','Audio','AV'};

ied_count_baseline = zeros(length(soz_types),length(mod_types));
ied_count_flicker = zeros(length(soz_types),length(mod_types));
ied_pct_change = zeros(length(soz_types),length(mod_types));

modal_est = zeros(length(soz_types),length(regions),length(mod_types));
modal_se = zeros(length(soz_types),length(regions),length(mod_types));
modal_baseline = zeros(length(soz_types),length(regions),length(mod_types));
modal_pvals = zeros(length(soz_types),length(regions),length(mod_types));
modal_ci_upper = zeros(length(soz_types),length(regions),length(mod_types));
modal_ci_lower = zeros(length(soz_types),length(regions),length(mod_types));
modal_df = zeros(length(soz_types),length(regions),length(mod_types));

%subs = unique(ied_table_wb.subjectID);

for kk0 = 1:length(soz_types)
    load([output_dir, '/IED_table_wb.mat'])
    pts_to_exc = unique(ied_table_wb.subjectID(~strcmp(ied_table_wb.soz_anat,soz_types{kk0})));
    ied_table_wb(contains(ied_table_wb.subjectID,pts_to_exc),:) = [];
    %ied_table_wb(strcmp(ied_table_wb.subjectID,subs{kk0}),:) = [];
    
    for kk = 1:length(regions)
        for kk2 = 1:length(mod_types)
            ied_table_temp = ied_table_wb(strcmp(ied_table_wb.modality,mod_types{kk2}),:);
            if kk == 1
                glme = fitglme(ied_table_temp, 'ied_count_mtl ~ 1 + baseline + (1|subjectID)',...
                    'Distribution','Poisson','Link','log');

            elseif kk == 2
                glme = fitglme(ied_table_temp, 'ied_count_pfc ~ 1 + baseline + (1|subjectID)',...
                    'Distribution','Poisson','Link','log');

            elseif kk==4
                glme = fitglme(ied_table_temp, 'ied_count_audio ~ 1 + baseline + (1|subjectID)',...
                    'Distribution','Poisson','Link','log');

                ied_means_audio(kk2,1) = mean(ied_table_temp.ied_count_audio(ied_table_temp.baseline == 0));
                ied_means_audio(kk2,2) = mean(ied_table_temp.ied_count_audio(ied_table_temp.baseline == 1));
            elseif kk==3
                glme = fitglme(ied_table_temp, 'ied_count_visual ~ 1 + baseline + (1|subjectID)',...
                    'Distribution','Poisson','Link','log');

                ied_means_visual(kk2,1) = mean(ied_table_temp.ied_count_visual(ied_table_temp.baseline == 0));
                ied_means_visual(kk2,2) = mean(ied_table_temp.ied_count_visual(ied_table_temp.baseline == 1));
            elseif kk==5
                glme = fitglme(ied_table_temp, 'ied_count_other ~ 1 + baseline + (1|subjectID)',...
                    'Distribution','Poisson','Link','log');
            end
            
            modal_pvals(kk0,kk,kk2) = glme.Coefficients.pValue(2);
            modal_est(kk0,kk,kk2) = glme.Coefficients.Estimate(2);
            modal_se(kk0,kk,kk2) = glme.Coefficients.SE(2);
            modal_ci_lower(kk0,kk,kk2) = glme.Coefficients.Lower(2);
            modal_ci_upper(kk0,kk,kk2) = glme.Coefficients.Upper(2);
            modal_baseline(kk0,kk,kk2) = exp(glme.Coefficients.Estimate(1));
            modal_df(kk0,kk,kk2) = glme.Coefficients.DF(2);
        end
    end
end

temp_struct = struct();
temp_struct.soz_types = soz_types;
temp_struct.modalities = mod_types;
temp_struct.regions = regions;
temp_struct.pvals = modal_pvals;
temp_struct.est = modal_est;
temp_struct.se = modal_se;
temp_struct.ci_upper = modal_ci_upper;
temp_struct.ci_lower = modal_ci_lower;
temp_struct.baseline = modal_baseline;
temp_struct.df = modal_df;

temp_struct.pct_est= -100*(exp(modal_est)-1);
temp_struct.pct_se = -100*(exp(modal_se)-1);
temp_struct.pct_ci_upper = -100*(exp(modal_ci_upper)-1);
temp_struct.pct_ci_lower = -100*(exp(modal_ci_lower)-1);

ied_results.region_modality_soztle = temp_struct;

%% Generate and save structure with per-session means (supplement A)
load([output_dir, '/ied_table_wb.mat'])
unique_sess = unique(ied_table_wb(:,[1,2,3]),'rows');

pct_change_all = nan(size(unique_sess,1),1);
base_all = nan(size(unique_sess,1),1);
fl_all = nan(size(unique_sess,1),1);

for kk = 1:size(unique_sess,1)
    tempinds = and(strcmp(ied_table_wb.subjectID,unique_sess.subjectID{kk}),strcmp(ied_table_wb.task_name,unique_sess.task_name{kk}));
    tempinds = and(tempinds, ied_table_wb.ses_nber == unique_sess.ses_nber(kk));
    subtable = ied_table_wb(tempinds,:);
    pct_change_all(kk) = 100*((nanmean(subtable.ied_count(subtable.baseline == 0))-nanmean(subtable.ied_count(subtable.baseline == 1)))./nanmean(subtable.ied_count(subtable.baseline == 1)));
    base_all(kk) = nanmean(subtable.ied_count(subtable.baseline == 1));
    fl_all(kk) = nanmean(subtable.ied_count(subtable.baseline == 0));
end

ied_by_sess.all_chan.pct = pct_change_all;
ied_by_sess.all_chan.baseline = base_all;
ied_by_sess.all_chan.flicker = fl_all;

save([output_dir, '/IED_by_session.mat'], 'ied_by_sess')

%% save data structure for stats
save([output_dir, '/IED_struct.mat'], 'ied_results')

disp('Results saved.')

function modalities = get_modalities(cond_arr)

    modalities = cell(size(cond_arr));
    
    for kk = 1:length(cond_arr)
        if contains(cond_arr{kk}, '-V')
            modalities{kk} = 'Visual';
        elseif contains(cond_arr{kk}, '-AV')
            modalities{kk} = 'AV';
        else
            modalities{kk} = 'Audio';
        end
    end
end
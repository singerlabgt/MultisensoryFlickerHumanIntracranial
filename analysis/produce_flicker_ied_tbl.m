%Creates table of IED counts along with all variables of interest, for the Poisson model.
%2024/02/26

function produce_flicker_ied_tbl(root_dir,output_dir)
    %% fetch metadata:
    metadata_tbl=readtable([root_dir '/FlickerStudyMetadata.xlsx'],'Sheet','Subjects','PreserveVariableNames',1);
    metadata_tbl=metadata_tbl(:,{'Subject_ID','Seizure focus broader classification'});

    %% fetch data for all sessions:
    %define sessions from fickerneuro and flickerfreq tasks:
    task_name='all';
    [~,all_sessions]=fetch_flicker_subjectIDs(root_dir,task_name);
    all_sessions(strcmp(all_sessions.task,'spep'),:)=[];

    %set some criteria:
    sample_rate=200; %sample rate of ied preprocessed data
    time_window=0.1; %time window to consider spikes to be part of same event
    max_channels=11; %maximum number of channels that can be involved in detecting given event

    %retrieve data:
    sessions=struct;
    for s=1:size(all_sessions,1)

        %define session information:
        sessions(s).subjectID=all_sessions{s,'sub'}{:};
        sessions(s).deidentified_subject_ID=metadata_tbl{strcmp(metadata_tbl.Subject_ID,sessions(s).subjectID),'Subject_ID'}{:};
        sessions(s).soz=metadata_tbl{strcmp(metadata_tbl.Subject_ID,sessions(s).subjectID),'Seizure focus broader classification'}{:};
        sessions(s).task=all_sessions{s,'task'}{:};
        sessions(s).ses=all_sessions{s,'ses'}{:};
        sessions(s).modality=all_sessions{s,'modality'}{:};

        %read ied data and clean:
        sessions(s).ied_tbl=readtable([root_dir '/stg-preproc/sub-' sessions(s).subjectID '/task-' sessions(s).task '/ses-' sessions(s).ses '/IED-preproc/sub-' sessions(s).subjectID '_allspikes.csv'],'Delimiter',','); %read detected ieds
        sessions(s).ied_tbl(logical(sessions(s).ied_tbl.predicted_class),:)=[]; %remove spikes considered not actual ieds
        sessions(s).ied_tbl(:,{'subject','clip_ids','clip','predicted_class'})=[]; %remove info we don't need
        sessions(s).ied_tbl.chan=sessions(s).ied_tbl.chan+1; %increment index of channels by 1 (0 in python is 1 in matlab)

        %convert spikes to events (i.e. a train of rapidly succeeding spikes, within time_window, is 1 event):
        %initialize event table:
        sessions(s).event_tbl=table('Size',[0 7],'VariableNames',{'start',...           %start of event, in time samples
                                                                  'start_mean',...      %mean start of event (i.e. mean of starts of all spikes constituting the event), in time samples
                                                                  'duration',...        %duration of the event (with last spike start + time_window), in time samples
                                                                  'spike_starts',...    %start times of each of the spikes in the event, in time samples
                                                                  'channels',...        %index of channels implicated in detecting each of the spikes in the event
                                                                  'num_spikes',...      %number of spikes in the event
                                                                  'num_channels'},...   %number of unique channels implicated in detecting spikes in this event
                                                                  'VariableTypes',{'double','double','double','cell','cell','double','double'});
        %fill-in event table:
        %STOPPED REVIEWING HERE...
        spike_list=1; %keeps track of index of spikes to include in given event; start with spike of index 1
        for sp=1:size(sessions(s).ied_tbl,1) %for each detected ied
            if sp==size(sessions(s).ied_tbl,1) || sessions(s).ied_tbl.start(sp+1)-sessions(s).ied_tbl.start(sp)>time_window*sample_rate %if we've reached end of ied tbl or the next spike is farther than time window, the current ied event is done and should be added to event table
                sessions(s).event_tbl(end+1,:)={sessions(s).ied_tbl.start(spike_list(1)),...
                                                mean(sessions(s).ied_tbl.start(spike_list)),...
                                                ((sessions(s).ied_tbl.start(spike_list(end))+time_window*sample_rate)-sessions(s).ied_tbl.start(spike_list(1))),...
                                                mat2cell(sessions(s).ied_tbl.start(spike_list)',1,length(spike_list)),...
                                                mat2cell(sessions(s).ied_tbl.chan(spike_list)',1,length(sessions(s).ied_tbl.chan(spike_list)')),...
                                                length(spike_list),...
                                                length(unique(sessions(s).ied_tbl.chan(spike_list)))};
                spike_list=sp+1; %start new list of index of spikes i.e. new event, starting with next spike
            elseif sessions(s).ied_tbl.start(sp+1)-sessions(s).ied_tbl.start(sp)<=time_window*sample_rate %if next spike is within time window
                spike_list=[spike_list,sp+1]; %include the next ied in the current event
            end
        end
        %sanity check that event table makes sense:

        sessions(s).event_tbl(sessions(s).event_tbl.num_channels>max_channels,:)=[]; %remove events that were detected by more than max_channels (considered noise)

        %get event start times:
        sessions(s).ied_event_times=sessions(s).event_tbl.start(:)'/sample_rate; %convert to absolute times in seconds

        %get information about trial times:
        sessions(s).trials=importdata([root_dir '/stg-preproc/sub-' sessions(s).subjectID '/task-' sessions(s).task '/ses-' sessions(s).ses '/sub-' sessions(s).subjectID '_stg-preproc_task-' sessions(s).task '_ses-' sessions(s).ses '_nat-beh.mat'],'trials');

        %get information about flicker modulation:
        sessions(s).ssep_amp_sig=readtable([root_dir '/stg-analyses/task-' sessions(s).task '/sub-' sessions(s).subjectID '/ses-' sessions(s).ses '/LFP/static_ent/LFP_pvalue_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
        sessions(s).ssep_amp_val=readtable([root_dir '/stg-analyses/task-' sessions(s).task '/sub-' sessions(s).subjectID '/ses-' sessions(s).ses '/LFP/static_ent/LFP_zscore_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);

        %get channel labels for events table:
        sessions(s).preproc_labels=readcell([root_dir '/stg-preproc/sub-' sessions(s).subjectID '/task-' sessions(s).task '/ses-' sessions(s).ses '/IED-preproc/sub-' sessions(s).subjectID '_eegdata_labels.csv']);
        sessions(s).preproc_labels=sessions(s).preproc_labels';
        
        sessions(s).sozchs={''};
        
        %ADD CODE TO CHECK THAT ALL PREPROC LABELS MATCH ALL MODULATION
        %ANALYSIS LABELS
        
        %fetch subject anat data:
        temp=fetch_subject_data(root_dir,{sessions(s).subjectID},'anat');
        sessions(s).anat=temp.anat;
    end
    
    %WRITE CODE TO CHECK FOR SUBCLINICAL SEIZURES

    %SHOW MEAN AND STANDARD DEVIATION OF SPIKES FOR STIM AND BASELINE, PER SUBJECT, COMBINING ALL CONDITIONS.

    %% define conditions, and regions of interest:
    %get list of conditions in order, for each of 2 tasks:
    flickerneuro_conditions={'5.5Hz-V','5.5Hz-AV','5.5Hz-A','40Hz-V','40Hz-AV','40Hz-A','80Hz-V','80Hz-AV','80Hz-A','R-V','R-AV','R-A'};
    if any(strcmp({sessions.task},'flickerfreq')) %if there are flickerfreq sessions
        temp=find(strcmp({sessions.task},'flickerfreq'));
        flickerfreq_conditions=sessions(temp(1)).ssep_amp_sig.Properties.VariableNames;
        flickerfreq_conditions=regexprep(flickerfreq_conditions,'-.','');
        [~,temp]=sort(arrayfun(@(x) str2num(x{:}),regexprep(flickerfreq_conditions,'Hz','')));
        flickerfreq_conditions=flickerfreq_conditions(temp);
        flickerfreq_conditions=[flickerfreq_conditions 'R']; %add random condition
    end

    %fetch freesurfer subregions for given regions:
    classifying_table=readtable('anat_cat.csv');
    regions=struct();
    regions.visual_regions=classifying_table.subregion_wmparc(strcmp(classifying_table.cog,'visual'));
    regions.audio_regions=classifying_table.subregion_wmparc(strcmp(classifying_table.cog,'auditory'));
    regions.MTL_regions=classifying_table.subregion_wmparc(strcmp(classifying_table.MTL_PFC,'MTL'));
    regions.PFC_regions=classifying_table.subregion_wmparc(strcmp(classifying_table.MTL_PFC,'PFC'));
   
    %% count number of ieds per trial:
    modch_percent_detectedIEDs=zeros(sum(strcmp({sessions.task},'flickerneuro')),length(flickerneuro_conditions)-3); %keep track of % IED events detected by group of modulated channels, for each session and condition.
    sozch_percent_detectedIEDs=zeros(1,sum(strcmp({sessions.task},'flickerneuro'))); %keep track of % IED events detected by group of soz channels, for each session and condition.
    regionch_percent_detectedIEDs=zeros(sum(strcmp({sessions.task},'flickerneuro')),length(fieldnames(regions)));
    flickerneuro_ses_nber=0;
    for s=1:length(sessions)
        
        %get overall spike rate across experiment:
        if isfield(sessions(s).trials.clinrecording,'SyncPulse')
            field_name='SyncPulse';
        else
            field_name='VisualOut';
        end
        sessions(s).overall_spikerate=size(sessions(s).event_tbl,1)/(length(sessions(s).trials.clinrecording.(field_name))/sessions(s).trials.clinrecording.sampleRate/60);

        %get average spike rate across baseline trials:
        baseline_trials_index=find(sessions(s).trials.trials_identities(:,1)==find(strcmp(sessions(s).trials.condition_code,'Baseline')));
        num_baseline_spikes=0;
        for tr=baseline_trials_index'
            num_baseline_spikes=num_baseline_spikes + sum(sessions(s).ied_event_times*sessions(s).trials.clinrecording.sampleRate>=sessions(s).trials.clinrecording.trials_timestamps(tr,1) & ...
                                 sessions(s).ied_event_times*sessions(s).trials.clinrecording.sampleRate<sessions(s).trials.clinrecording.trials_timestamps(tr,2)); %get number of ieds that occurred within that baseline trial;
        end
        sessions(s).baseline_ied_rate=(num_baseline_spikes/(length(baseline_trials_index)*10))*60;

        if strcmp(sessions(s).task,'flickerfreq')
            conditions=flickerfreq_conditions;

            %calculate proportion of spikes that occur in stim (vs surrounding 5x2s
            %baseline) trials:
            sessions(s).proportion_stim_spikes=nan(10,length(conditions));
            total_trials=sum(~contains(sessions(s).trials.condition_code(sessions(s).trials.trials_identities(:,1)),'occluded')); %get total number of trials that are not occluded trials
            for cond=1:length(conditions)
                trials_index=find(sessions(s).trials.trials_identities(:,1)==find(strcmp(sessions(s).trials.condition_code,[conditions{cond} '-' sessions(s).modality]))); %find index of trials for that condition
                for tr=1:length(trials_index)
                    %calculate proportion of spikes that occurred during the stim
                    %part of trial:
                    temp_stimieds=sum(sessions(s).ied_event_times*sessions(s).trials.clinrecording.sampleRate>=sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr),1) & ...
                                     sessions(s).ied_event_times*sessions(s).trials.clinrecording.sampleRate<sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr),2)); %number of ieds occurring during stim trial

                    if ismember(trials_index(tr),[1 2]) %if we're within first 2 trials of experiment
                        temp_baseline_trials=1:5;
                    elseif ismember(trials_index(tr),total_trials-1:total_trials) %if we're within last 2 trials of experiment
                        temp_baseline_trials=total_trials-4:total_trials;
                    else
                        temp_baseline_trials=trials_index(tr)-2:trials_index(tr)+2;
                    end

                    temp_baselineieds=0;
                    for baseline_tr=1:length(temp_baseline_trials)
                        temp_baselineieds=temp_baselineieds+sum(sessions(s).ied_event_times*sessions(s).trials.clinrecording.sampleRate>=sessions(s).trials.clinrecording.trials_timestamps(temp_baseline_trials(baseline_tr),2) & ...
                                         sessions(s).ied_event_times*sessions(s).trials.clinrecording.sampleRate<(sessions(s).trials.clinrecording.trials_timestamps(temp_baseline_trials(baseline_tr),2)+2*sessions(s).trials.clinrecording.sampleRate));
                    end

                    sessions(s).proportion_stim_spikes(tr,cond)=temp_stimieds/(temp_stimieds+temp_baselineieds); %proportion of stim ieds (adjusted for longer stim trial) over whole length of stim+baseline trial.
                end
            end

            %count spikes occurring for each channel, condition (stim or corresponding baseline), and trial:
            sessions(s).spikes_ch_tr=nan(size(sessions(s).ssep_amp_sig,1),length(conditions),10,2); %4D matrix for ch, condition, whether stim or not, and trial number
            total_trials=sum(~contains(sessions(s).trials.condition_code(sessions(s).trials.trials_identities(:,1)),'occluded')); %get total number of trials that are not occluded trials
            for cond=1:length(conditions)
                for ch=1:size(sessions(s).spikes_ch_tr,1)
                    ied_ch=find(strcmp(sessions(s).preproc_labels,regexprep(sessions(s).ssep_amp_sig.Properties.RowNames{ch},'-.+','')));
                    ied_event_times=sessions(s).event_tbl;
                    temp=arrayfun(@(x) any(x{:}(1)==ied_ch),ied_event_times.channels,'UniformOutput',false);
                    ied_event_times=sessions(s).ied_event_times([temp{:}]);

                    trials_index=find(sessions(s).trials.trials_identities(:,1)==find(strcmp(sessions(s).trials.condition_code,[conditions{cond} '-' sessions(s).modality]))); %find index of trials for that condition
                    for tr=1:length(trials_index)
                        sessions(s).spikes_ch_tr(ch,cond,tr,1)=sum(ied_event_times*sessions(s).trials.clinrecording.sampleRate>=sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr),1) & ...
                            ied_event_times*sessions(s).trials.clinrecording.sampleRate<sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr),2));

                        if ismember(trials_index(tr),[1 2]) %if we're within first 2 trials of experiment
                            temp_baseline_trials=1:5;
                        elseif ismember(trials_index(tr),total_trials-1:total_trials) %if we're within last 2 trials of experiment
                            temp_baseline_trials=total_trials-4:total_trials;
                        else
                            temp_baseline_trials=trials_index(tr)-2:trials_index(tr)+2;
                        end

                        temp_baselineieds=0;
                        for baseline_tr=1:length(temp_baseline_trials)
                            temp_baselineieds=temp_baselineieds+sum(ied_event_times*sessions(s).trials.clinrecording.sampleRate>=sessions(s).trials.clinrecording.trials_timestamps(temp_baseline_trials(baseline_tr),2) & ...
                                ied_event_times*sessions(s).trials.clinrecording.sampleRate<(sessions(s).trials.clinrecording.trials_timestamps(temp_baseline_trials(baseline_tr),2)+2*sessions(s).trials.clinrecording.sampleRate));
                        end


                        sessions(s).spikes_ch_tr(ch,cond,tr,2)=temp_baselineieds;
                    end
                end
            end

        elseif strcmp(sessions(s).task,'flickerneuro')
            conditions=flickerneuro_conditions;
            flickerneuro_ses_nber=flickerneuro_ses_nber+1;

            %calculate proportion of spikes that occur in stim (vs subsequent
            %baseline) trials:
            sessions(s).proportion_stim_spikes=nan(15,length(conditions));
            for cond=1:length(conditions)
                trials_index=find(sessions(s).trials.trials_identities(:,1)==find(strcmp(sessions(s).trials.condition_code,conditions{cond}))); %find index of trials for that condition
                for tr=1:length(trials_index)
                    %calculate proportion of spikes that occurred during the stim
                    %part of trial:
                    sessions(s).proportion_stim_spikes(tr,cond)=sum(sessions(s).ied_event_times*sessions(s).trials.clinrecording.sampleRate>=sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr),1) & ...
                        sessions(s).ied_event_times*sessions(s).trials.clinrecording.sampleRate<sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr),2))/... %number of spikes occurring in the stim part of trial
                        sum(sessions(s).ied_event_times*sessions(s).trials.clinrecording.sampleRate>=sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr),1) & ...
                        sessions(s).ied_event_times*sessions(s).trials.clinrecording.sampleRate<sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr)+1,2)); %number of spikes occurring in the stim+baseline part of trial
                end
            end

            %count spikes occurring for each channel, condition (stim or corresponding baseline), and trial:
            sessions(s).spikes_ch_tr=nan(size(sessions(s).ssep_amp_sig,1),length(conditions),15,2); %4D matrix for ch, condition, whether stim or not, and trial number
            for cond=1:length(conditions)
                for ch=1:size(sessions(s).spikes_ch_tr,1)
                    ied_ch=find(strcmp(sessions(s).preproc_labels,regexprep(sessions(s).ssep_amp_sig.Properties.RowNames{ch},'-.+','')));
                    ied_event_times=sessions(s).event_tbl;
                    temp=arrayfun(@(x) any(x{:}(1)==ied_ch),ied_event_times.channels,'UniformOutput',false);
                    ied_event_times=sessions(s).ied_event_times([temp{:}]);

                    trials_index=find(sessions(s).trials.trials_identities(:,1)==find(strcmp(sessions(s).trials.condition_code,conditions{cond}))); %find index of trials for that condition
                    for tr=1:length(trials_index)
                        sessions(s).spikes_ch_tr(ch,cond,tr,1)=sum(ied_event_times*sessions(s).trials.clinrecording.sampleRate>=sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr),1) & ...
                            ied_event_times*sessions(s).trials.clinrecording.sampleRate<sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr),2));
                        sessions(s).spikes_ch_tr(ch,cond,tr,2)=sum(ied_event_times*sessions(s).trials.clinrecording.sampleRate>=sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr)+1,1) & ...
                            ied_event_times*sessions(s).trials.clinrecording.sampleRate<sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr)+1,2));
                    end
                end
            end
        end
    end

    modch_percent_detectedIEDs=array2table(modch_percent_detectedIEDs,'VariableNames',flickerneuro_conditions(1:end-3));
    sozch_percent_detectedIEDs=sozch_percent_detectedIEDs';
    regionch_percent_detectedIEDs=array2table(regionch_percent_detectedIEDs,'VariableNames',fieldnames(regions)');

    %% organize data into one table- for analysis with Poisson Generalized Linear Mixed effects Model (GLMM):

    data_values=zeros([10000000 7]);
    data_labels=cell([10000000 8]);
    row_nber=0;
    for s=1:length(sessions)
        if strcmp(sessions(s).task,'flickerfreq')
            conditions=strcat(flickerfreq_conditions,'-',sessions(s).modality);
        elseif strcmp(sessions(s).task,'flickerneuro')
            conditions=flickerneuro_conditions;
        end

        s_init=row_nber;
        total_trials=sum(~contains(sessions(s).trials.condition_code(sessions(s).trials.trials_identities(:,1)),'occluded')); %get total number of trials that are not occluded trials
        for ch=1:size(sessions(s).ssep_amp_sig,1)
            ch_init=row_nber;
            ied_ch=find(strcmp(sessions(s).preproc_labels,regexprep(sessions(s).ssep_amp_sig.Properties.RowNames{ch},'-.+','')));
            ied_event_times=sessions(s).event_tbl;
            temp=arrayfun(@(x) any(x{:}(1)==ied_ch),ied_event_times.channels,'UniformOutput',false);
            ied_event_times=sessions(s).ied_event_times([temp{:}]);
            for cond=1:length(conditions)
                cond_init=row_nber;
                trials_index=find(sessions(s).trials.trials_identities(:,1)==find(strcmp(sessions(s).trials.condition_code,conditions{cond}))); %find index of trials for that condition
                for tr=1:length(trials_index)
                        for baseline=[0 1]
                            row_nber=row_nber+1;
                            data_values(row_nber,4:5)=[baseline,sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr),1)/sessions(s).trials.clinrecording.sampleRate];
                            if ~baseline
                                data_values(row_nber,6)=sum(ied_event_times*sessions(s).trials.clinrecording.sampleRate>=sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr),1) & ...
                                                            ied_event_times*sessions(s).trials.clinrecording.sampleRate<sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr),2));
                            elseif baseline
                                if strcmp(sessions(s).task,'flickerneuro')
                                    data_values(row_nber,6)=sum(ied_event_times*sessions(s).trials.clinrecording.sampleRate>=sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr)+1,1) & ...
                                                            ied_event_times*sessions(s).trials.clinrecording.sampleRate<sessions(s).trials.clinrecording.trials_timestamps(trials_index(tr)+1,2));
                                elseif strcmp(sessions(s).task,'flickerfreq')
                                    if ismember(trials_index(tr),[1 2]) %if we're within first 2 trials of experiment
                                        temp_baseline_trials=1:5;
                                    elseif ismember(trials_index(tr),total_trials-1:total_trials) %if we're within last 2 trials of experiment
                                        temp_baseline_trials=total_trials-4:total_trials;
                                    else
                                        temp_baseline_trials=trials_index(tr)-2:trials_index(tr)+2;
                                    end

                                    temp_baselineieds=0;
                                    for baseline_tr=1:length(temp_baseline_trials)
                                        temp_baselineieds=temp_baselineieds+sum(ied_event_times*sessions(s).trials.clinrecording.sampleRate>=sessions(s).trials.clinrecording.trials_timestamps(temp_baseline_trials(baseline_tr),2) & ...
                                            ied_event_times*sessions(s).trials.clinrecording.sampleRate<(sessions(s).trials.clinrecording.trials_timestamps(temp_baseline_trials(baseline_tr),2)+2*sessions(s).trials.clinrecording.sampleRate));
                                    end

                                    data_values(row_nber,6)=temp_baselineieds;
                                end
                            end
                        end
                end
                if ~contains(conditions{cond},'R-')
                    data_values(cond_init+1:cond_init+row_nber-cond_init,7:8)=repmat([sessions(s).ssep_amp_sig{ch,conditions{cond}},sessions(s).ssep_amp_val{ch,conditions{cond}}],row_nber-cond_init,1);
                else
                    data_values(cond_init+1:cond_init+row_nber-cond_init,7:8)=repmat([NaN,NaN],row_nber-cond_init,1);
                end
                data_labels(cond_init+1:cond_init+row_nber-cond_init,8)=repmat({conditions{cond}},row_nber-cond_init,1);
            end
            ied_ch=find(strcmp(sessions(s).preproc_labels,regexprep(sessions(s).ssep_amp_sig.Properties.RowNames{ch},'-.+','')));
            ied_event_times=sessions(s).event_tbl;
            temp=arrayfun(@(x) any(x{:}(1)==ied_ch),ied_event_times.channels,'UniformOutput',false);
            ied_event_times=sessions(s).ied_event_times([temp{:}]);

            data_values(ch_init+1:ch_init+row_nber-ch_init,2:3)=repmat([ismember(regexprep(sessions(s).ssep_amp_sig.Properties.RowNames{ch},'-.+',''),sessions(s).sozchs)...
                                                                        length(ied_event_times)/size(sessions(s).event_tbl,1)],row_nber-ch_init,1);
            data_labels(ch_init+1:ch_init+row_nber-ch_init,5:6)=repmat({regexprep(sessions(s).ssep_amp_sig.Properties.RowNames{ch},'-.+','')...
                                                                        sessions(s).anat.electrodes_info.anatlabels{strcmp(sessions(s).anat.electrodes_info.labels,regexprep(sessions(s).ssep_amp_sig.Properties.RowNames{ch},'-.+','')),'fs_aparcaseg'}},row_nber-ch_init,1);
            for temp_field={'visual_regions','audio_regions','MTL_regions','PFC_regions'}
                if contains(data_labels{ch_init+1,6},regions.(temp_field{:}))
                    data_labels(ch_init+1:ch_init+row_nber-ch_init,7)=repmat({temp_field{:}},row_nber-ch_init,1);
                    break;
                end
            end
            if isempty(data_labels{ch_init+1,7})
                data_labels(ch_init+1:ch_init+row_nber-ch_init,7)=repmat({'other'},row_nber-ch_init,1);
            end
        end
        data_values(s_init+1:s_init+row_nber-s_init,1)=repmat(sessions(s).overall_spikerate,row_nber-s_init,1);
        data_labels(s_init+1:s_init+row_nber-s_init,1:4)=repmat({sessions(s).deidentified_subject_ID,sessions(s).task,sessions(s).ses,sessions(s).soz},row_nber-s_init,1);
    end

    data_tbl=[cell2table(data_labels(logical(data_values(:,5)),:),'VariableNames',{'subjectID','task_name','ses_nber','soz_anat','ch','ch_anat','ch_region','cond'}) array2table(data_values(logical(data_values(:,5)),:),'VariableNames',{'session_ied_rate','in_soz','ch_ied_proportion','baseline','trial_start_time','ied_count','mod_sig','mod_amp'})];

    %combine some anatomy:
    temp=strings(size(data_tbl,1),1);
    temp(strcmp(data_tbl.ch_anat,'') | contains(data_tbl.ch_anat,{'Ventricle','choroid'}) | endsWith(data_tbl.ch_anat,'Vent'))='other';
    temp(~(strcmp(data_tbl.ch_anat,'') | contains(data_tbl.ch_anat,{'Ventricle','choroid'}) | endsWith(data_tbl.ch_anat,'Vent')))=regexprep(data_tbl{~(strcmp(data_tbl.ch_anat,'') | contains(data_tbl.ch_anat,{'Ventricle','choroid'}) | endsWith(data_tbl.ch_anat,'Vent')),'ch_anat'},{'-rh','-lh','Right-','Left-'},'');
    temp=array2table(temp,'VariableNames',{'ch_anat_combined'});
    data_tbl=[data_tbl(:,1:6),temp, data_tbl(:,7:end)];

    %writetable(data_tbl,[figures_dir '/IED_count_table_all-sessions.csv']);
    data_tbl(:,{'ch_anat','ch_anat_combined','session_ied_rate','in_soz','ch_ied_proportion'})=[];
    writetable(data_tbl,[output_dir '/IED_count_table_all-sessions.csv']);
    
end

%Evaluate flicker LFP phase locking value significance and amplitude.
%2024/02/25

function evaluate_ssep_plv(fnames)
    
    %create parent output dir if does not exist:
    if ~exist([fnames.analysis_folder,'/LFP/static_ent'],'dir')
        mkdir([fnames.analysis_folder,'/LFP/static_ent']);
    end
    
    for ref_method={'Laplacian'}
        %retrieve data:
        trials=importdata([fnames.preprocdata_folder '/sub-' fnames.subjectID '_stg-preproc_task-' fnames.task '_ses-' fnames.ses '_nat-beh.mat'],'trials'); %get trials data
        data=importdata([fnames.preproc_LFPdata '/sub-' fnames.subjectID '_stg-preproc_task-' fnames.task '_ses-' fnames.ses '_nat-trialsdata-ref' ref_method{:} '-preproc.mat'],'data_ref_preproc'); %get preprocessed/segmented LFP data
        
        %define conditions of interest:
        conditions_of_interest=trials.condition_code(contains(trials.condition_code,'Hz') & ~contains(trials.condition_code,'occluded') & ~contains(trials.condition_code,'min'));
        conditions_of_interest=conditions_of_interest(ismember(conditions_of_interest,trials.condition_code(trials.trials_identities(:,1))));
        conditions_of_interest=sort(conditions_of_interest);
        
        t=data.time{1};
        t_index=find(t>=0 & t<=10);
        fsample=trials.clinrecording.sampleRate;
        
        pvalue_table=zeros(length(data.label),length(conditions_of_interest));
        plv_table=zeros(length(data.label),length(conditions_of_interest));
        for cond=1:size(plv_table,2) %for each condition
            stim_freq=str2num(regexprep(conditions_of_interest{cond},'Hz.+','')); %obtain stim frequency
            flicker_signal=sin(2*pi*stim_freq*t); %simulate flicker signal
       
            for ch=1:size(plv_table,1) %for each channel
                angle_diff_t=[];
                for tr=find(data.trialinfo==find(strcmp(trials.condition_code,conditions_of_interest(cond))))' %for each trial of this condition
                    tmp_neuralsignal=data.trial{tr}(ch,:); %get neural signal
                    [s,f,time,c] = xspectrogram(tmp_neuralsignal(t_index),flicker_signal(t_index),hanning(fsample/2),0,[stim_freq-1.5 stim_freq stim_freq+1.5],fsample);
                    angle_diff_t=[angle_diff_t;exp(1i*angle(c(2,:)))]; %get phase angle difference between flicker signal and neural signal, for trial from 0-10s
                end
                plv_t=mean(angle_diff_t); %take the angle difference per time bin, averaged across trials
                plv_table(ch,cond)=abs(mean(plv_t)); %get average of average phase angle, across time
                pvalue_table(ch,cond)=circ_rtest(angle(plv_t)); %calculate significance value, using circular statistics
            end
        end
    
        %save results:
        plv_table=array2table(plv_table,'RowNames',data.label,'VariableNames',conditions_of_interest); %make matrix into table
        pvalue_table=array2table(pvalue_table,'RowNames',data.label,'VariableNames',conditions_of_interest); %make matrix into table
        
        writetable(plv_table,[fnames.analysis_folder,'/LFP/static_ent/LFP_plv-amp_table_ref' ref_method{:} '.csv'],'WriteRowNames',1);
        writetable(pvalue_table,[fnames.analysis_folder,'/LFP/static_ent/LFP_plv-pval_table_ref' ref_method{:} '.csv'],'WriteRowNames',1);
    end
end

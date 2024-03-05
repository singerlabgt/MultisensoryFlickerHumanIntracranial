%This function creates simulated 40Hz flicker data based on the SPEP
%response and the superposition hypothesis.
%NEED TO ADJUST DATA VARIABLE SO THAT ITS FIELDS MATCH THE KINDS OF FIELDS
%WE GET FOR FLICKER FIELDTRIP DATA.
%2024/02/25

function superposition_modeling_data(fnames)
    
    %fetch data:
    load(fnames.trials_filename); %fetch trials data
    load([fnames.preproc_LFPdata '/sub-' fnames.subjectID '_stg-analysis_task-' fnames.task '_ses-' fnames.ses '_nat-spepERP-refLaplacian.mat'],'ERP_results_ref_preproc'); %fetch preprocessed data
    
    %create flicker trials structure:
    flicker_trials=struct();
    flicker_trials.condition_code=["40Hz-V";"40Hz-AV";"40Hz-A";"Baseline"];
    flicker_trials.clinrecording.sampleRate=trials.clinrecording.sampleRate;
    
    %create simulated data structure (in fieldtrip format):
    data=struct('fsample',[],'trialinfo',[],'trial',[],'time',[],'label',[],'ref','Laplacian');
    data.fsample=flicker_trials.clinrecording.sampleRate;
    data.trialinfo=[];
    for i=1:length(flicker_trials.condition_code)
        data.trialinfo=[data.trialinfo;repmat(i,15,1)];
    end
    data.time{1,1}=ERP_results_ref_preproc{1}.time(1):1/data.fsample:11-ERP_results_ref_preproc{1}.time(1);
    data.time{1,1}=data.time{1,1}(1:end-1);
    for i=1:length(data.trialinfo)
        data.time{1,i}=data.time{1,1};
    end
    data.label=ERP_results_ref_preproc{1}.label;
    
    %simulate the data for each trial:
    sample_mult_factor=40/gcd(40,data.fsample);
    for i=1:length(data.trialinfo) %for each trial
        flicker_condition=flicker_trials.condition_code(data.trialinfo(i)); %determine which flicker condition we're simulating
        if strcmp(flicker_condition,'Baseline')
            spep_condition='occluded_AV';
        else
            spep_condition=strsplit(flicker_condition,'-');
            spep_condition=spep_condition(end);
        end
        spep_condition_nber=0;
        for temp=1:4
            if ERP_results_ref_preproc{temp}.trialinfo==find(strcmp(trials.condition_code,spep_condition))
                spep_condition_nber=temp;
            end
        end
        
        for j=1:length(data.label) %for each channel
            current_timeseries=zeros(1,length(data.time{1,1})*sample_mult_factor); %initialize trial LFP at value 0, upsampled
            current_index=1;
            %index_tracking=current_index;
            for k=1:400 %simulate 10s 40Hz flicker trial (40x10=400 pulses)
                current_spep=ERP_results_ref_preproc{spep_condition_nber}.trial(randsample(size(ERP_results_ref_preproc{spep_condition_nber}.trial,1),1),j,:); %randomly pick 1 SPEP from that condition and channel
                current_spep=reshape(current_spep,[1 size(current_spep,3)]);
                current_spep=interp(current_spep,sample_mult_factor); %resample spep to same sample rate as current_timeseries
                current_timeseries(current_index:current_index+length(current_spep)-1)=current_timeseries(current_index:current_index+length(current_spep)-1)+current_spep;
                current_index=current_index+(1/40)*data.fsample*sample_mult_factor;
                %index_tracking=[index_tracking current_index];
            end
            current_timeseries=downsample(current_timeseries,sample_mult_factor);
            data.trial{1,i}(j,:)=current_timeseries;
        end
    end
    
    %create output parent folder if does not exist:
    if ~exist(fnames.analysis_folder,'dir')
        mkdir(fnames.analysis_folder);
    end
    
    %save simulated data:
    trials=flicker_trials;
    save([fnames.analysis_folder '/sub-' fnames.subjectID '_stg-analysis_model-superposition_ses-' fnames.ses '_nat-trialsdata-refLaplacian.mat'],'trials','data','-v7.3');
end

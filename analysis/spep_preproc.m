%Segments and preprocessed LFP from spep task.
%2024/02/25

function spep_preproc(fnames)
    %load data:
    load(fnames.trials_filename);

    %organize trial data into fieldtrip format:
    analysis_trials=[trials.clinrecording.trials_timestamps(1,1) trials.clinrecording.trials_timestamps(end,end)];
    analysis_trials(:,1)=analysis_trials(:,1)-60*trials.clinrecording.sampleRate; %take 0.2s before each trial
    analysis_trials(:,2)=analysis_trials(:,2)+60*trials.clinrecording.sampleRate;
    analysis_trials(:,3)=-60*trials.clinrecording.sampleRate; %specify when trial actually starts
    analysis_trials(:,4)=0;

    %preprocess recordings from 1s trials:
    %fetch raw data:
    cfg=[];
    cfg.datafile=fnames.LFP_filename;
    cfg.hdr= ft_read_header(fnames.LFP_filename);
    channels=fetch_channels(fnames,'-chNoisy'); %fetch channels we're interested in analyzing
    cfg.channel=channels.channel; %retrive data from channels we're interested in analyzing
    cfg.trl=analysis_trials;
    raw_data=ft_preprocessing(cfg);
    raw_data.label=channels.label; %correct channel labels if needed
    
    %save this raw segmented data:
    outputFolder=[fnames.preprocdata_folder '/LFP'];
    if ~exist([outputFolder '/spep_erp'],'dir')
        mkdir([outputFolder '/spep_erp']);
    end
    
    %apply preprocessing:
    %define parameters to apply to re-referenced LFP segments:
    cfg1=[];
    cfg1.hpfilter       = 'yes'; %highpass filter %FIGURE OUT WHY HIGH PASS FILTER DOESN'T WORK
    cfg1.hpfreq         = 0.1; %to get rid of DC shift
    cfg1.hpfilttype='but';
    cfg1.hpfiltord=4;
    
    %format trials information into fieldtrip format:
    analysis_trials=trials.clinrecording.trials_timestamps;
    analysis_trials(:,1)=uint32(analysis_trials(:,1)-0.25*trials.clinrecording.sampleRate); %take 0.2s before each trial %NEED TO DO UINT32 HERE?
    analysis_trials(:,2)=uint32(analysis_trials(:,2)+0.25*trials.clinrecording.sampleRate);
    analysis_trials(:,3)=-0.25*trials.clinrecording.sampleRate; %specify when trial actually starts
    analysis_trials(:,4)=trials.trials_identities(:,1);
    cfg2=[];
    cfg2.trl=analysis_trials;
    
    %to apply baseline correction
    cfg3.demean='yes';
    cfg3.baselinewindow=[-Inf 0]; %baseline correction will be applied using mean of signal 0.25s before start of pulse DOES USING -INF MEAN GO TO THE SMALLEST TIME IN TRIAL?
    
    %apply re-referencing and preprocessing:
    for ref_method={'Laplacian'}
        switch ref_method{:}
            case 'bipolar' % apply bipolar referencing:
                disp('Applying bipolar montage...');
                bipolar=create_bipolar_montage(raw_data.label);
                data_ref=ft_apply_montage(raw_data,bipolar);
            case 'Laplacian'
                disp('Applying Laplacian montage...');
                Laplacian=create_Laplacian_montage(raw_data.label);
                data_ref=ft_apply_montage(raw_data,Laplacian);
        end
        
        disp('Preprocessing data...');
        data_ref_preproc_pre = ft_preprocessing(cfg1,data_ref);
        clear data_ref;
                
        disp('Resegmenting data...');
        data_ref_preproc_reseg=ft_redefinetrial(cfg2,data_ref_preproc_pre);
        clear data_ref_preproc_pre;
        
        disp('Applying baseline correction...');
        data_ref_preproc=ft_preprocessing(cfg3,data_ref_preproc_reseg);
        clear data_ref_preproc_reseg;
        
        %save preprocessed data:
        disp('Saving preprocessed data...');
        data_ref_preproc.ref=ref_method;
        save([outputFolder '/spep_erp/sub-' fnames.subjectID '_stg-preproc_task-' fnames.task '_ses-' fnames.ses '_nat-trialsdata-ref' ref_method{:} '-preproc.mat'],'data_ref_preproc','-v7.3'); %save this data
        clear data_ref_preproc;
    end
end



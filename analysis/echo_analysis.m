%Performs analysis on the persistence of oscillatory activity following offset of flicker stimulus.
%Requirements:
%MATLAB version 2021b or later because of the pyrunfile function.
%Python. To make sure have proper python version/environment, type in: pyversion(<path to proper python.exe file>).
%2024/02/25

function echo_analysis(fnames)
    %load trials data, and preprocessed/semented data:
    data=importdata([fnames.preproc_LFPdata 'sub-' fnames.subjectID '_stg-preproc_task-' fnames.task '_ses-' fnames.ses '_nat-trialsdata-refLaplacian-preproc.mat']);
    trials=importdata(fnames.trials_filename);
    
    %clean line noise:
    cfg=[];
    cfg.dftfilter='yes';
    cfg.dftfreq=60;
    cfg.dftreplace='neighbour';
    data=ft_preprocessing(cfg,data);
    
    %define conditions of interest:
    if strcmp(fnames.task,'flickerneuro') %for the flickerneuro task
        conditions={'5.5Hz-V','5.5Hz-AV','5.5Hz-A','40Hz-V','40Hz-AV','40Hz-A','80Hz-V','80Hz-AV','80Hz-A'};
    elseif strcmp(fnames.task,'flickerfreq') %for the flickerfreq task
        conditions=unique(trials.condition_code(data.trialinfo));
        conditions(contains(conditions,{'occluded','Baseline','R'}))=[];
        [~,temp_index]=sort(arrayfun(@(x) str2num(x),regexprep(conditions,'Hz-.+','')));
        conditions=cellstr(conditions(temp_index))';
    end
    
    %run analysis:
    echo_results=cell(1,length(conditions)); %initialize array of results
    for c=conditions %for each condition of interest
        %keep track of conditions, channel labels and sample rate:
        echo_results{strcmp(conditions,c)}.condition=c;
        echo_results{strcmp(conditions,c)}.label=data.label;
        echo_results{strcmp(conditions,c)}.sampleRate=trials.clinrecording.sampleRate;
        
        %get 10s flicker ERP response, averaged across trials:
        cfg=[];
        cfg.trials=find(data.trialinfo==find(strcmp(trials.condition_code,c)));
        cfg.keeptrials='no';
        temp=ft_timelockanalysis(cfg,data);
        echo_results{strcmp(conditions,c)}.evoked=temp.avg;
    
        %run echo analysis on evoked responses, using script adapted from Lerousseau et al. J. Neurosci. 2021:
        [echo_results{strcmp(conditions,c)}.evoked_filtered,echo_results{strcmp(conditions,c)}.evoked_hil,echo_results{strcmp(conditions,c)}.evoked_hil_z,echo_results{strcmp(conditions,c)}.thresh,...
            echo_results{strcmp(conditions,c)}.evoked_hil_rec,echo_results{strcmp(conditions,c)}.evoked_hil_rec_bin,echo_results{strcmp(conditions,c)}.threshold,echo_results{strcmp(conditions,c)}.good_channels,...
            echo_results{strcmp(conditions,c)}.n_cycles,echo_results{strcmp(conditions,c)}.onsets]...
            =pyrunfile('echo_analysis.py',... %python script to run
            ["EVOKED_FILTERED","EVOKED_HIL","EVOKED_HIL_Z","thresh","EVOKED_HIL_REC","EVOKED_HIL_REC_BIN","threshold","GOOD_CHANNELS","N_CYCLES","ONSETS"],... %variables from python script to return, at end of script
            sfreq=echo_results{strcmp(conditions,c)}.sampleRate,trial_duration=10,prestim_duration=1,poststim_duration=1,condition_freq=str2double(regexprep(c,'Hz-.+','')),EVOKED=echo_results{strcmp(conditions,c)}.evoked); %varibles to provide to python script
        
        %get rid of data we do not need:
        for fieldname={'evoked_filtered','evoked_hil','evoked_hil_z','thresh','evoked_hil_rec','evoked_hil_rec_bin','threshold','good_channels','n_cycles','onsets'}
            echo_results{strcmp(conditions,c)}.(fieldname{:})=double(echo_results{strcmp(conditions,c)}.(fieldname{:}));
        end
    end
    
    %save results:
    save([fnames.analysis_folder '/LFP/static_ent/sub-' fnames.subjectID '_stg-analysis_task-' fnames.task '_ses-' fnames.ses '_nat-echo-analysis-preproc_cleaned-linenoise.mat'],'echo_results','-v7.3');
end



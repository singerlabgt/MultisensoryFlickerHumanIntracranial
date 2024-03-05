%Performs PSD analysis for all trials in a flicker session, and saves data.
%2024/02/25

function PSD_results=run_PSD_flicker(data,trials)
    %initialize PSD_results variable:
    PSD_results.label=data.label;
    PSD_results.condition=repmat({''},1,length(trials.condition_code));
    
    % Reformat data into chronux format (samples*trials):
    disp('Preparing data for Chronux...');
    data_for_chronux=cell(length(PSD_results.label),length(PSD_results.condition));
    [~,temp]=min(abs(data.time{1}-10)); %to deal with cases where trial does not have an index for time 10s (because of sample rate not being an integer); CORRECT WAY TO DO THIS?
    time_index=find(data.time{1}==0):temp; %take only the trial itself, i.e. not the extra seconds before and after
    for i=1:length(trials.condition_code) %for each condition
        trials_of_interest=find(data.trialinfo==i); %find trials of interest
        if ~isempty(trials_of_interest)
            disp('Preparing data from 1 condition...');
            PSD_results.condition{i}=char(trials.condition_code(i));
            for k=trials_of_interest' %for each trial
                for j=1:length(PSD_results.label) %for each electrode
                    data_for_chronux{j,i}(:,end+1)=data.trial{k}(j,time_index);
                end
            end
        end
    end

    %perform power spectral density analysis per condition (using Chronux toolbox):
    %set parameters:
    disp('Performing PSD analysis...');
    params.Fs = trials.clinrecording.sampleRate;
    params.fpass = [2 100];
    params.err = [1 0.05];
    params.trialave = 0;
    params.tapers=[3 5];
    %perform analysis:
    PSD_results.data=cell(size(data_for_chronux));
    %note: mtspectrumc uses dpss from MATLAB (not from FieldTrip, otherwise error)
    for i=1:length(PSD_results.condition) %for each condition
        if ~isempty(PSD_results.condition{i})
            disp(['Processing 1 condition out of ' num2str(sum(cellfun(@(x) ~isempty(x),PSD_results.condition))) ' conditions']);
            for j=1:length(PSD_results.label) %for each electrode
                [S,f,Serr]=mtspectrumc(data_for_chronux{j,i},params);
                PSD_results.data{j,i}{1}=S';
                PSD_results.data{j,i}{2}=Serr;
                PSD_results.data{j,i}{3}=f;
            end
        end
    end
end

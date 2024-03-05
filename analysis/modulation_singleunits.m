%Estimates flicker modulation of single units, using PSTH plots an circular statistics.
%2024/02/25

function modulation_singleunits(fnames)

    %CONTINUE WORKING ON THIS...
    %SHOW VISUAL AND AUDIO RESPONSE FOR SAME UNIT, OR SHOW A SUBSET OF
    %CONDITIONS FOR SAME UNIT, SO READER CAN COMPARE. MOREOVER, RANDOM SHOULD
    %BE ALIGNED TO PULSE ONSET.

    %load trials data and unit data:
    trials=importdata([fnames.preprocdata_folder '/sub-' fnames.subjectID '_stg-preproc_task-' fnames.task '_ses-' fnames.ses '_nat-beh.mat'],'trials'); %contains trial data;
    unit=importdata([fnames.preprocdata_folder '/single-unit-preproc/sub-' fnames.subjectID '_stg-preproc_task-' fnames.task '_ses-' fnames.ses '_nat-unit-data.mat'],'unit'); %contains single unit data

    %define output folder:
    outputFolder=[fnames.analysis_folder '/single-unit'];
    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end

    %organize trials data into 2 repeats of pulse:
    [conditionIndex, data]=converttimestamps3(trials,'BRrecording',2,1);

    %only keep a random set of 15 baseline trials (so number of baseline trials matches number of stim trials):
    for condition=find(contains(conditionIndex,'Baseline'))
        random_set_of_baseline_trials=datasample(unique(data{condition}(:,4)),15,'Replace',false);
        data{condition}=data{condition}(ismember(data{condition}(:,4),random_set_of_baseline_trials),:);
    end

    %convert unit structure:
    spike_original=struct();
    spike_original.label={unit.name};
    spike_original.timestamp=arrayfun(@(x) x.spikeTimes*30,unit,'UniformOutput',0); %converted timestamps back to 1/30000 (Combinato gives results in ms).
    spike_original.waveform=arrayfun(@(x) {reshape(x.spikeWaveforms,[1,64,size(x.spikeWaveforms,2)])},unit); %NEED TO CHECK WHETHER COMBINATO OUTPUT WAVEFORMS AT ORIGINAL SAMPLE RATE
    
    %convert to FieldTrip format:
    spike_trials=data;
    cfg=[];
    cfg.timestampspersecond=30000;
    for i=1:length(data)
        spike_trials{i}(:,3:end)=[];
        spike_trials{i}(:,1)=(spike_trials{i}(:,1)-1)*(30000/trials.BRrecording.sampleRate)+1;
        spike_trials{i}(:,2)=spike_trials{i}(:,2)*(30000/trials.BRrecording.sampleRate);
        spike_trials{i}(:,3)=0;
        cfg.trl=spike_trials{i};
        spike_trials{i}=ft_spike_maketrials(cfg,spike_original);
    end

    %calculate PSTH:
    psth_data=spike_trials;
    cfg=[];
    cfg.outputunit='rate';
    cfg.keeptrials='no';
    for i=1:length(psth_data)
        cfg.binsize=psth_data{i}.trialtime(1,2)/20;
        psth_data{i}=ft_spike_psth(cfg,psth_data{i});
    end
    
    %plot 1 plot per unit, with all conditions plotted in same figure
    for u=1:length(unit) %for each unit
        figure('units','normalized','outerposition',[0 0 1 1]);
        tiledlayout(3,4,'TileSpacing','none','Padding','none');
        tile_order=[2 3 4 6 7 8 10 11 12];
        condition={'5.5Hz-V','5.5Hz-AV','5.5Hz-A','40Hz-V','40Hz-AV','40Hz-A','80Hz-V','80Hz-AV','80Hz-A'}; %define flickerneuro conditions
        nexttile(1);
        plot(unit(u).spikeWaveforms,'Color',[0 0 0 0.05]);
        title(['Unit name: ' unit(u).name '; Nber spikes: ' num2str(unit(u).numberSpikes)],'Interpreter','none');
        for i=1:length(condition) %for each condition
            nexttile(tile_order(i));
            unit_label=unit(u).name;
            plot_flicker_psth(psth_data,conditionIndex,condition{i},unit_label,NaN);
        end

        %save unit plots:
        print([outputFolder '/sub-' fnames.subjectID '_ses-' fnames.ses '_singleUnitModulation.ps'],'-dpsc','-bestfit','-append');
        close all;
    end
end

%2024/02/26

%% define directories:
root_dir=define_flicker_root_dir;
figures_dir=[root_dir '/stg-analyses/NatureComm2024-figures'];

%% fetch single unit data:
%CONTINUE WORKING ON THIS...
%RANDOM SHOULD BE ALIGNED TO PULSE ONSET?

metadata_tbl=readtable([root_dir '/FlickerStudyMetadata.xlsx'],'Sheet','Sessions','PreserveVariableNames',1);
metadata_tbl=metadata_tbl(metadata_tbl.("Single units included in NatureComm2024?")==1,:);

%specify anatomical location of microwires:
micro_anat=table2array(metadata_tbl(:,{'Subject_ID','abbreviated location of single units'}));
micro_anat=[micro_anat(:,1) repmat({1:8},size(micro_anat,1),1) micro_anat(:,2)];
to_remove=[];
to_add=[];
for i=1:size(micro_anat,1)
    if contains(micro_anat(i,3),',') %means there are multiple single unit locations (each location has 8 microwires)
        to_remove=[to_remove, i];
        temp_locations=strsplit(micro_anat{i,3},', ');
        to_add=[repmat(micro_anat(i,1),length(temp_locations),1) num2cell(reshape(1:length(temp_locations)*8,8,2)',2) temp_locations'];
    end
end
micro_anat(to_remove,:)=[];
micro_anat=[micro_anat;to_add];
micro_anat=sortrows(micro_anat,1);

%load trials data and unit data for all sessions:
unit_data=struct();
unit_data_nber=0;
for sub_nber=1:size(metadata_tbl)
    unit_data_nber=unit_data_nber+1;
    unit_data(unit_data_nber).subjectID=metadata_tbl{sub_nber,'Subject_ID'}{:};
    unit_data(unit_data_nber).task=metadata_tbl{sub_nber,'Experiment'}{:};
    unit_data(unit_data_nber).ses=['0' num2str(metadata_tbl{sub_nber,'Session'})];
    unit_data(unit_data_nber).trials=importdata([root_dir '/stg-preproc/sub-' unit_data(unit_data_nber).subjectID '/task-flickerneuro/ses-' unit_data(unit_data_nber).ses '/sub-' unit_data(unit_data_nber).subjectID '_stg-preproc_task-flickerNeuro_ses-' unit_data(unit_data_nber).ses '_nat-beh.mat'],'trials'); %contains trial data;
    unit_data(unit_data_nber).unit=importdata([root_dir '/stg-preproc/sub-' unit_data(unit_data_nber).subjectID '/task-flickerneuro/ses-' unit_data(unit_data_nber).ses '/single-unit-preproc/sub-' unit_data(unit_data_nber).subjectID '_stg-preproc_task-flickerneuro_ses-' unit_data(unit_data_nber).ses '_nat-unit-data.mat'],'unit'); %contains single unit data
    
    %assign unit anatomical location:
    temp_subject=unit_data(unit_data_nber).subjectID;
    for u=1:length(unit_data(unit_data_nber).unit)
        temp_anat=micro_anat(strcmp(micro_anat(:,1),temp_subject),:);
        unit_data(unit_data_nber).unit(u).anat=temp_anat{arrayfun(@(x) ismember(unit_data(unit_data_nber).unit(u).ch,x{:}),temp_anat(:,2)),3};
    end
end

%display number of isolated units:
su_mu=zeros(1,2);
for u_d_nber=1:length(unit_data)
    su_mu(1)=su_mu(1)+sum(strcmp({unit_data(u_d_nber).unit.status},'SU'));
    su_mu(2)=su_mu(2)+sum(strcmp({unit_data(u_d_nber).unit.status},'MU'));
end
disp(['Total of ' num2str(sum(su_mu)) ' units isolated: ' num2str(su_mu(1)) ' single units, ' num2str(su_mu(2)) ' multi-units.']);

%show histogram of spike rates:
spike_rates=[];
for i=1:length(unit_data)
    for j=1:length(unit_data(i).unit)
        spike_rates=[spike_rates unit_data(i).unit(j).numberSpikes/(length(unit_data(i).trials.BRrecording.SyncPulse)/unit_data(i).trials.BRrecording.sampleRate/60)];
    end
end

histogram(spike_rates,50)
title('Distribution of spike rates (spikes/min)');


%% transform single unit data and align to pulses:

for i=1:length(unit_data)
    
    %convert unit structure to FieldTrip format:
    unit_data(i).spike_original=struct();
    unit_data(i).spike_original.label={unit_data(i).unit.name};
    unit_data(i).spike_original.timestamp=arrayfun(@(x) x.spikeTimes*30,unit_data(i).unit,'UniformOutput',0); %converted timestamps back to 1/30000 (Combinato gives results in ms).
    unit_data(i).spike_original.waveform=arrayfun(@(x) {reshape(x.spikeWaveforms,[1,64,size(x.spikeWaveforms,2)])},unit_data(i).unit); %NEED TO CHECK WHETHER COMBINATO OUTPUT WAVEFORMS AT ORIGINAL SAMPLE RATE
    
    %organize trials data into 2 repeats of pulse:
    [unit_data(i).trials_converted.condition_index,unit_data(i).trials_converted.data]=converttimestamps3(unit_data(i).trials,'BRrecording',2,1);

    %only keep a random set of 15 baseline trials (so number of baseline trials matches number of stim trials):
    for condition=find(contains(unit_data(i).trials_converted.condition_index,'Baseline'))
        random_set_of_baseline_trials=datasample(unique(unit_data(i).trials_converted.data{condition}(:,4)),15,'Replace',false);
        unit_data(i).trials_converted.data{condition}=unit_data(i).trials_converted.data{condition}(ismember(unit_data(i).trials_converted.data{condition}(:,4),random_set_of_baseline_trials),:);
    end
    
    unit_data(i).spike_trials=unit_data(i).trials_converted.data;
    cfg=[];
    cfg.timestampspersecond=30000;
    for j=1:length(unit_data(i).trials_converted.data)
        unit_data(i).spike_trials{j}(:,3:end)=[];
        unit_data(i).spike_trials{j}(:,1)=(unit_data(i).spike_trials{j}(:,1)-1)*(30000/unit_data(i).trials.BRrecording.sampleRate)+1;
        unit_data(i).spike_trials{j}(:,2)=unit_data(i).spike_trials{j}(:,2)*(30000/unit_data(i).trials.BRrecording.sampleRate);
        unit_data(i).spike_trials{j}(:,3)=0;
        cfg.trl=unit_data(i).spike_trials{j};
        unit_data(i).spike_trials{j}=ft_spike_maketrials(cfg,unit_data(i).spike_original);
    end
    
    unit_data(i).psth_data=unit_data(i).spike_trials;
    cfg=[];
    cfg.outputunit='rate';
    cfg.keeptrials='no';
    for j=1:length(unit_data(i).psth_data)
        cfg.binsize=unit_data(i).psth_data{j}.trialtime(1,2)/20;
        unit_data(i).psth_data{j}=ft_spike_psth(cfg,unit_data(i).psth_data{j});
    end
    
    unit_data(i).spike_trials_radians=unit_data(i).spike_trials;
    for j=1:length(unit_data(i).spike_trials_radians)
        freq=strsplit(unit_data(i).trials_converted.condition_index{j},'-');
        freq=regexprep(freq{1},'Hz','');
        freq=str2double(freq);
        for k=1:length(unit_data(i).spike_trials_radians{j}.time)
            unit_data(i).spike_trials_radians{j}.time{k}=unit_data(i).spike_trials_radians{j}.time{k}*2*pi*freq;
        end
    end
    
end


%% exclude single units that have too few spikes:
conditions={'5.5Hz-V','5.5Hz-AV','5.5Hz-A','40Hz-V','40Hz-AV','40Hz-A','80Hz-V','80Hz-AV','80Hz-A'};

row_nber=0;
units_to_exclude=[];
for ses=1:length(unit_data)
    for u=1:length(unit_data(ses).unit)
        row_nber=row_nber+1;
        %determine if all conditions have more than 20% of the bins with
        %null spike rate- eliminate those units based on this criterium:
        temp=[];
        for cond=1:length(conditions)
            temp=[temp sum(~unit_data(ses).psth_data{strcmp(unit_data(ses).trials_converted.condition_index,conditions{cond})}.avg(u,:))/size(unit_data(ses).psth_data{strcmp(unit_data(ses).trials_converted.condition_index,conditions{cond})}.avg,2)];
        end
        if all(temp>0.2)
            units_to_exclude=[units_to_exclude row_nber];
        end
    end
end

disp(['Eliminating ' num2str(length(units_to_exclude)) ' units based on selection criteria; left with ' num2str(sum(arrayfun(@(x) length(x{:}),{unit_data.unit}))-length(units_to_exclude)) ' units.']);

%spell out how many SU and MU, and where (hippo or cingulate) we included
%in analysis:
hippoCing_SuMu=zeros(2,2); %matrix keep track of number of SU vs MU (columns) in hippo vs cingulate (rows)
unit_nber=0;
for unit_data_nber=1:length(unit_data)
    for u=1:length(unit_data(unit_data_nber).unit)
        unit_nber=unit_nber+1;
        if ~ismember(unit_nber,units_to_exclude)
            if strcmp(unit_data(unit_data_nber).unit(u).anat,'hippocampus') && strcmp(unit_data(unit_data_nber).unit(u).status,'SU')
                hippoCing_SuMu(1,1)=hippoCing_SuMu(1,1)+1;
            elseif strcmp(unit_data(unit_data_nber).unit(u).anat,'cingulate') && strcmp(unit_data(unit_data_nber).unit(u).status,'SU')
                hippoCing_SuMu(2,1)=hippoCing_SuMu(2,1)+1;
            elseif strcmp(unit_data(unit_data_nber).unit(u).anat,'hippocampus') && strcmp(unit_data(unit_data_nber).unit(u).status,'MU')
                hippoCing_SuMu(1,2)=hippoCing_SuMu(1,2)+1;
            elseif strcmp(unit_data(unit_data_nber).unit(u).anat,'cingulate') && strcmp(unit_data(unit_data_nber).unit(u).status,'MU')
                hippoCing_SuMu(2,2)=hippoCing_SuMu(2,2)+1;
            end
        end
    end
end

hippoCing_SuMu

%% draw PSTHs for example single units and conditions:

%define examples:
examples=readtable([figures_dir '/NatureComm2024_examples.xlsx']);
examples=table2array(examples(startsWith(examples.example,'Figure_S5'),:));
[~,temp_index]=sort(examples(:,1));
examples=examples(temp_index,:);
examples=examples(:,[2 5 6 1]);

for u=1:size(examples,1)
    figure_making('width',6.5,'height',1.5);
    tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
    
    ses_nber=find(strcmp({unit_data.subjectID},examples{u,1}));
    unit_nber=find(strcmp({unit_data(ses_nber).unit.name},examples{u,2}));
    condition_freq=examples{u,3};
    
    %plot all spikes:
    nexttile;
    plot(unit_data(ses_nber).unit(unit_nber).spikeWaveforms,'Color',[0 0 0 0.01]);
    hold on;
    plot(mean(unit_data(ses_nber).unit(unit_nber).spikeWaveforms'),'k');
    title([num2str(round(unit_data(ses_nber).unit(unit_nber).numberSpikes/(length(unit_data(ses_nber).trials.BRrecording.SyncPulse)/unit_data(ses_nber).trials.BRrecording.sampleRate/60),0)) ' spikes/min'],'FontSize',7);
    box off;
    axis tight;
    set(gca,'Visible','off');
    ax=gca;
    add_scale_bars(ax,0,20,0.01);
    ax.Title.Visible='on';
    
    %plot psth for each modality at that frequency:
    conditions=strcat([condition_freq '-'],{'V','AV','A'});
    max_spikerate=0;
    for cond=1:length(conditions)
        temp=max(unit_data(ses_nber).psth_data{strcmp(unit_data(ses_nber).trials_converted.condition_index,conditions{cond})}.avg(unit_nber,:));
        if temp>max_spikerate
            max_spikerate=temp;
        end
    end
    max_spikerate=max_spikerate+max_spikerate/5;
    for cond=conditions
        nexttile;
        plot_flicker_psth(unit_data(ses_nber).psth_data,unit_data(ses_nber).trials_converted.condition_index,cond{:},examples{u,2},max_spikerate);
        title([]);
        ax=gca;
        if contains(cond,'5.5Hz')
            ax.XTick=round(ax.XTick([1 3])*1000)/1000;
        else
            ax.XTick=round(ax.XTick([1 3])*1000,1)/1000;
        end
        ax.XTickLabel=string(ax.XTick*1000);
        if ~contains(cond,'AV')
            xlabel([]);
        end
        if ~contains(cond,'-V')
            ylabel([]);
        else
            ylabel('Spike rate (/s)');
        end
        temp=gca;
        
        condition_VS=circ_r(unit_data(ses_nber).spike_trials_radians{strcmp(unit_data(ses_nber).trials_converted.condition_index,cond)}.time{unit_nber}');
        text(0+range(temp.XLim)/20,temp.YLim(end)-range(temp.YLim)/10,['VS: ' num2str(round(condition_VS,2))]);
        condition_RS=2*length(unit_data(ses_nber).spike_trials_radians{strcmp(unit_data(ses_nber).trials_converted.condition_index,cond)}.time{unit_nber})*condition_VS^2;
        text(temp.XLim(end)-range(temp.XLim)/2.5,temp.YLim(end)-range(temp.YLim)/10,['RS: ' num2str(round(condition_RS,2))]);
        
        temp1=strsplit(cond{:},'-');
        comparator=[temp1{1} '-R' temp1{2}];
        random_VS=circ_r(unit_data(ses_nber).spike_trials_radians{strcmp(unit_data(ses_nber).trials_converted.condition_index,comparator)}.time{unit_nber}');
        text(0+range(temp.XLim)/20,temp.YLim(1)+range(temp.YLim)/15,['VS: ' num2str(round(random_VS,2))]);
        random_RS=2*length(unit_data(ses_nber).spike_trials_radians{strcmp(unit_data(ses_nber).trials_converted.condition_index,comparator)}.time{unit_nber})*random_VS^2;
        text(temp.XLim(end)-range(temp.XLim)/2.5,temp.YLim(1)+range(temp.YLim)/15,['RS: ' num2str(round(random_RS,2))]);
    end
    
    figure_making('operation','save','filename',[figures_dir '/' examples{u,4} '.pdf']);
end

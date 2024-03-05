%Shows examples of raw LFP, segmented LFP and PSD for given ch and condition.
%Includes parts of panels in figures 1C, 1E, 2A, 2B, S2A, S2B, S3A.
%2024/02/26

%% set directory paths:
root_dir=define_flicker_root_dir;
figures_dir=[root_dir '/stg-analyses/NatureComm2024-figures'];

%% define details of examples, fetch data and preprocess it:

%define which examples we'd like to show:
examples=readtable([figures_dir '/NatureComm2024_examples.xlsx']);
examples=table2array(examples(startsWith(examples.example,{'Figure_1','Figure_2','Figure_S2','Figure_S3'}),:));
examples=[examples(:,2:end) examples(:,1)];

%determine which sessions we need to retrieve data from:
sessions_info=examples(:,1:3);
sessions_info=unique(join(string(sessions_info),',')); %find unique set of sessions
sessions_info=cellfun(@(x) strsplit(x,','),sessions_info,'UniformOutput',false);
sessions_info=vertcat(sessions_info{:});

%fetch data and  minimize space it takes:
sessions=struct();
for s=1:size(sessions_info,1)
    disp(['Fetching data for session ' num2str(s) ' out of ' num2str(size(sessions_info,1))]);
    sessions(s).subjectID=sessions_info{s,1};
    sessions(s).task=sessions_info{s,2};
    sessions(s).ses=sessions_info{s,3};
    sessions(s).trials=importdata([root_dir '/stg-preproc/sub-' sessions(s).subjectID '/task-' sessions(s).task '/ses-' sessions(s).ses '/sub-' sessions(s).subjectID '_stg-preproc_task-' sessions(s).task '_ses-' sessions(s).ses '_nat-beh.mat'],'trials');
    sessions(s).electrodes_info=importdata([root_dir '/stg-preproc/sub-' sessions(s).subjectID '/anat/electrode-models/sub-' sessions(s).subjectID '_electrodes_info.mat'],'electrodes_info');
    
    %data below takes space, so only keep what we need:
    ex_nbers=find(strcmp(examples(:,1),sessions(s).subjectID) & strcmp(examples(:,2),sessions(s).task) & strcmp(examples(:,3),sessions(s).ses));
    temp=importdata([root_dir '/stg-preproc/sub-' sessions(s).subjectID '/task-' sessions(s).task '/ses-' sessions(s).ses '/LFP/static_ent/sub-' sessions(s).subjectID '_stg-analysis_task-' sessions(s).task '_ses-' sessions(s).ses '_nat-psd-refLaplacian.mat'],'PSD_results_ref_preproc');
    sessions(s).clinLFP_contacts=extract_clinLFP_labels(temp.label);
    
    temp_labels=regexprep(examples(ex_nbers,4),'\d+$','');
    temp_labels=vertcat(sessions(s).clinLFP_contacts(ismember({sessions(s).clinLFP_contacts.depth_electrode_name},temp_labels)).channel_names);
    
    chs_del=~ismember(temp.label,temp_labels);
    conds_del=~ismember(temp.condition,[examples{ex_nbers,5},{'Baseline'}]);
    temp.label(chs_del)=[];
    temp.condition(conds_del)=[];
    temp.data(chs_del,:)=[];
    temp.data(:,conds_del)=[];
    sessions(s).psd_data=temp;
    clear temp;

    temp=importdata([root_dir '/stg-preproc/sub-' sessions(s).subjectID '/task-' sessions(s).task '/ses-' sessions(s).ses '/LFP/static_ent/sub-' sessions(s).subjectID '_stg-preproc_task-' sessions(s).task '_ses-' sessions(s).ses '_nat-trialsdata-refLaplacian-preproc.mat'],'data_ref_preproc');    
    chs_del=~startsWith(temp.label,strcat(examples(ex_nbers,4),'-'));
    tr_del=~ismember(temp.trialinfo,find(ismember(sessions(s).trials.condition_code,[examples{ex_nbers,5},{'Baseline'}])));
    temp.sampleinfo(tr_del,:)=[];
    temp.trialinfo(tr_del)=[];
    temp.trial(tr_del)=[];
    for tr=1:length(temp.trial)
        temp.trial{tr}(chs_del,:)=[];
    end
    temp.time(tr_del)=[];
    temp.label(chs_del)=[];
    sessions(s).lfp_data=temp;
    clear temp;
end

%preprocess data:
examples_data=struct();
for ex=1:size(examples,1)

    %define example info:
    examples_data(ex).subjectID=examples{ex,1};
    examples_data(ex).task=examples{ex,2};
    examples_data(ex).ses=examples{ex,3};
    examples_data(ex).condition=examples{ex,5};
    examples_data(ex).type=examples{ex,6};

    %find which session we are going to use data from:
    s=find(strcmp(sessions_info(:,1),examples{ex,1}) & strcmp(sessions_info(:,2),examples{ex,2}) & strcmp(sessions_info(:,3),examples{ex,3}));

    examples_data(ex).channel=sessions(s).psd_data.label(startsWith(sessions(s).psd_data.label,[examples{ex,4} '-']));
    
    %pick example LFP trace and clean ground noise:
    temp=find(sessions(s).lfp_data.trialinfo==find(strcmp(sessions(s).trials.condition_code,examples_data(ex).condition)));
    temp=temp(6);
    cfg=[];
    cfg.channel=examples_data(ex).channel;
    cfg.bsfilter='yes';
    cfg.bsfreq=[58 62];
    cfg.trials=temp;
    temp=ft_preprocessing(cfg,sessions(s).lfp_data);
    examples_data(ex).lfp_example_cond=temp;
    
    %show evoked potential averaged across trials:
    cfg=[];
    cfg.channel=examples_data(ex).channel;
    cfg.trials=find(sessions(s).lfp_data.trialinfo==find(strcmp(sessions(s).trials.condition_code,examples_data(ex).condition)));
    cfg.keeptrials='no';
    temp=ft_timelockanalysis(cfg,sessions(s).lfp_data);
    examples_data(ex).lfp_average_cond=temp;

    %average LFP over 2 cycles of stim for condition and baseline:
    %NEED TO CHECK:
    [temp_index, temp_data]=converttimestamps3(sessions(s).trials,'clinrecording',2,1); %convert timestamps to timestamps for 2 cycles %NEED TO WORK ON THIS
    control_condition=strsplit(examples_data(ex).condition,'-');
    control_condition=[control_condition{1} '-Baseline'];
    analysis_trials=[];
    for cond=[{examples_data(ex).condition} control_condition]
        temp=temp_data{strcmp(temp_index,cond)};
        temp(:,3)=0;
        temp(:,4)=find(strcmp(temp_index,cond));
        analysis_trials(end+1:end+size(temp,1),:)=temp;
    end
    analysis_trials=uint32(analysis_trials);
    
    data_withTrials={};
    element_nber=0;
    for i=1:length(sessions(s).lfp_data.trial)
        cfg=[];
        cfg.trials=i;
        cfg.trl=analysis_trials(analysis_trials(:,1)>=sessions(s).lfp_data.sampleinfo(i,1) & analysis_trials(:,2)<=sessions(s).lfp_data.sampleinfo(i,2),:);
        if ~isempty(cfg.trl)
            element_nber=element_nber+1;
            data_withTrials{element_nber}=ft_redefinetrial(cfg,sessions(s).lfp_data);
        end
    end
    cfg=[];
    data_withTrials_sum=ft_appenddata(cfg,data_withTrials{:});
    
    cfg=[];
    cfg.channel=examples_data(ex).channel;
    cfg.trials=find(data_withTrials_sum.trialinfo==find(strcmp(temp_index,examples_data(ex).condition)));
    cfg.keeptrials='yes';
    examples_data(ex).lfp_psth_cond=ft_timelockanalysis(cfg,data_withTrials_sum);

    cfg=[];
    cfg.channel=examples_data(ex).channel;
    cfg.trials=find(data_withTrials_sum.trialinfo==find(strcmp(temp_index,control_condition)));
    cfg.keeptrials='yes';
    examples_data(ex).lfp_psth_control=ft_timelockanalysis(cfg,data_withTrials_sum);

end


%% draw plots:
%FOR FIGURE_S3A_3_EP, ERROR PATCH DOES NOT SHOW

for ex=1:length(examples_data)
    
    %find which session we are going to use data from:
    s=find(strcmp(sessions_info(:,1),examples_data(ex).subjectID) & strcmp(sessions_info(:,2),examples_data(ex).task) & strcmp(sessions_info(:,3),examples_data(ex).ses));

    %define frequency, modality, and corresponding color:
    temp=strsplit(examples_data(ex).condition,'-');
    freq_interest=str2num(regexprep(temp{1},'Hz',''));
    mod_interest=temp{2};
    line_color=condition_color(mod_interest);
    
    %plot example LFP trace (or averaged LFP trace):
    if startsWith(examples(ex,7),{'Figure_1','Figure_2','Figure_S2'}) 
        figure_making('width',3.2,'height',0.25);
        plot(examples_data(ex).lfp_example_cond.time{1},examples_data(ex).lfp_example_cond.trial{1},'Color','k','LineWidth',1);
        %plot(examples_data(ex).lfp_average_cond.time,examples_data(ex).lfp_average_cond.avg,'Color','k','LineWidth',1);
        xlim([-0.1 0.4]);
        temp=gca;
        temp=[min(temp.Children.YData(temp.Children.XData>=-0.1 & temp.Children.XData<=0.5)) max(temp.Children.YData(temp.Children.XData>=-0.1 & temp.Children.XData<=0.5))];
        plot_flicker(freq_interest,temp,10,line_color,0);
        ylim([temp(1) temp(2)+((temp(2)-temp(1))/5)]);
        box off;
        xline(0,'--','LineWidth',1,'Label','','LabelHorizontalAlignment','left');
        set(gca,'visible','off');
        
        ax=gca;
        ax.Position=[0.1 0.01 0.89 0.98];
        
        a2=annotation('line',[ax.Position(1) ax.Position(1)+0.1*(ax.Position(3)/range(ax.XLim))],[ax.Position(2) ax.Position(2)],'Linewidth',1);
        
        temp=num2str(range(ax.YLim));
        temp2=regexp(temp,'\.');
        temp=floor(range(ax.YLim)/(10^(temp2-2)))*(10^(temp2-2));
        a2=annotation('line',[ax.Position(1) ax.Position(1)],[ax.Position(2) ax.Position(2)+temp*(ax.Position(4)/range(ax.YLim))],'Linewidth',1);
        
        figure_making('operation','save','filename',[figures_dir '/' examples{ex,7} '_LFP.pdf']);
    end

    %plot LFP PSTH:
    figure_making('width',1,'height',1);
    erp_result=reshape(examples_data(ex).lfp_psth_cond.trial,size(examples_data(ex).lfp_psth_cond.trial,1),size(examples_data(ex).lfp_psth_cond.trial,3));
    plot(examples_data(ex).lfp_psth_cond.time,mean(erp_result),'Color',line_color,'LineWidth',0.5);
    hold on;
    x=[examples_data(ex).lfp_psth_cond.time(1:end-1) examples_data(ex).lfp_psth_cond.time(end-1:-1:1)];
    temp_patchvalues=[mean(erp_result(:,1:end-1))-std(erp_result(:,1:end-1))/sqrt(size(erp_result,1)) mean(erp_result(:,end-1:-1:1))+std(erp_result(:,end-1:-1:1))/sqrt(size(erp_result,1))];
    if ismember(ex,[3,4,8,9,10,12]) %need to do this because we have 2 NaN values in middle of temp (and so patch does fill):
        x(length(x)/2:length(x)/2+1)=[];
        temp_patchvalues(length(temp_patchvalues)/2:length(temp_patchvalues)/2+1)=[];
    end
    p=patch(x,temp_patchvalues,line_color,'FaceAlpha',0.2,'EdgeColor','none');
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    ax=gca;
    plot_flicker(freq_interest,ax.YLim,1/freq_interest*2,line_color,0,10);
    ylabel('Potential (\muV)');
    xlabel('Time (ms)');
    
    erp_result=reshape(examples_data(ex).lfp_psth_control.trial,size(examples_data(ex).lfp_psth_control.trial,1),size(examples_data(ex).lfp_psth_control.trial,3));
    plot(examples_data(ex).lfp_psth_control.time,mean(erp_result),'Color','k','LineWidth',0.5);
    x=[examples_data(ex).lfp_psth_control.time(1:end-1) examples_data(ex).lfp_psth_control.time(end-1:-1:1)];
    p=patch(x,[mean(erp_result(:,1:end-1))-std(erp_result(:,1:end-1))/sqrt(size(erp_result,1)) mean(erp_result(:,end-1:-1:1))+std(erp_result(:,end-1:-1:1))/sqrt(size(erp_result,1))],'k','FaceAlpha',0.2,'EdgeColor','none');
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    
    axis tight;
    ax=gca;
    ax.XTickLabel=arrayfun(@(x) num2str(str2num(x{:})*1000),ax.XTickLabel,'UniformOutput',false);
    box off;
    
    figure_making('operation','save','filename',[figures_dir '/' examples{ex,7} '_EP.pdf']);
    
    %plot PSD:
    figure_making('width',1,'height',1);
    plot_PSD('Baseline',examples_data(ex).channel,sessions(s).psd_data,sessions(s).psd_data.label,sessions(s).psd_data.condition,'k',1,0,1); %plot control condition for contact of interest
    hold on;
    plot_PSD(examples_data(ex).condition,examples_data(ex).channel,sessions(s).psd_data,sessions(s).psd_data.label,sessions(s).psd_data.condition,line_color,1,0,1); %plot condition
    %title(['PSD (' num2str(freq_interest) 'Hz vs R)'],'FontWeight','normal');
    
    switch freq_interest
        case 5.5
            xlim([0.5 10.5])
        case 40
            xlim([35 45]);
        case 80
            xlim([75 85]);
    end
    l=xline(freq_interest,'--','LineWidth',1);
    l.Alpha=0.5;
    p=gca;
    temp_ylim=p.YLim;
    %legend(p.Children([7,4]),{'R-AV',condition});
    ax=gca;
    %ax.XAxis.FontSize=9;
    xlabel('Frequency (Hz)');
    %ax.YAxis.FontSize=9;
    ylabel('Power (log_1_0(dB))');
    box off;
    
    %add legend:
    text(freq_interest+0.5,ax.YLim(end)-range(ax.YLim)/6,['\color[rgb]{' num2str(line_color) '}' examples_data(ex).condition '' newline '\color[rgb]{0 0 0}Baseline'],'FontSize',7);
    
    figure_making('operation','save','filename',[figures_dir '/' examples{ex,7} '_PSD.pdf']);
    
    %plot individual channels PSDs:
    if startsWith(examples(ex,7),'Figure_2')
        figure_making('width',5,'height',0.5);
        temp=strsplit(examples_data(ex).channel{:},'-');
        temp=regexprep(temp{1},'\d+$','');
        depth_electrode_contacts=sessions(s).clinLFP_contacts(strcmp({sessions(s).clinLFP_contacts.depth_electrode_name},temp)).channel_names;
        tiledlayout(1,15,'Padding','none','TileSpacing','compact');
        max_y_value=0;
        
        for i=1:length(depth_electrode_contacts)
            temp=strsplit(depth_electrode_contacts{i},'-');
            temp=regexp(temp{1},'\d','match');
            if length(temp)>1
                temp={strcat(temp{:})};
            end
            temp=str2num(temp{:});
            if ex==8
                nexttile(16-temp); %so order is flipped
            else
                nexttile(temp);
            end
            set(gca,'fontsize',0.1);
            plot_PSD(examples_data(ex).condition,depth_electrode_contacts{i},sessions(s).psd_data,sessions(s).psd_data.label,sessions(s).psd_data.condition,line_color,1,0,1); %plot condition
            ax = gca;
            set(ax,'XTick',[], 'YTick', []);
            box off;
            set(gca,'Visible','off');
            xlim([freq_interest-2 freq_interest+2]);
            temp=gca;
            median_value=median(temp.Children(3).YData(temp.Children(3).XData>=freq_interest-2 & temp.Children(3).XData<=freq_interest+2));
            ylim([median_value-0.4 median_value+0.4]);
            l=xline(freq_interest,'--','LineWidth',1);
            l.Alpha=0.5;
            temp=max(temp.Children(3).YData(temp.Children(3).XData>=freq_interest-2 & temp.Children(3).XData<=freq_interest+2))-median_value;
            if temp>max_y_value
                max_y_value=temp;
            end
        end
        
        for i=1:length(depth_electrode_contacts)
            temp=strsplit(depth_electrode_contacts{i},'-');
            temp=regexp(temp{1},'\d','match');
            if length(temp)>1
                temp={strcat(temp{:})};
            end
            temp=str2num(temp{:});
            if ex==2
                nexttile(16-temp); %so order is flipped
            else
                nexttile(temp);
            end
            temp=gca;
            ylim([temp.YLim(1)+(temp.YLim(2)-temp.YLim(1))/2-max_y_value temp.YLim(1)+(temp.YLim(2)-temp.YLim(1))/2+max_y_value]);
        end
        
        figure_making('operation','save','filename',[figures_dir '/' examples{ex,7} '_PSD_contacts.pdf']);
    end
end


%% fetch imaging and display contats of interest:
% _ load postop CT, and preop T1 transformed to postop CT space, and
% electrode model for postop CT space.
% _ both imaging in "smooth display" mode
% _ postop CT on top of preop MRI, with 15% opacity
% _ both imaging rotated by given degrees on given axes (so that "head is
% straight), and can see a few of the depth electrode contacts



%% retrieve all data:

%fetch all flicker sessions:
root_dir=define_flicker_root_dir;
figures_dir=[root_dir '/stg-analyses/NatureComm2024-figures'];
[~,all_sessions]=fetch_flicker_subjectIDs(root_dir,'all');
all_sessions(strcmp(all_sessions.task,'spep'),:)=[];

selected_examples=readtable([figures_dir '/NatureComm2024_examples.xlsx']);
selected_examples=table2array(selected_examples(startsWith(selected_examples.example,'Figure_S8'),:));
[~,index]=sort(selected_examples(:,1));
selected_examples(index,:);
selected_examples=selected_examples(:,[2 3 4 6 5 1]);

all_sessions=all_sessions(ismember(join(string(table2array(all_sessions(:,[1 2 3]))),';'),join(string(selected_examples(:,[1 2 3])),';')),:);


%% find highest SSEP in data:

flickerneuro_conditions={'5.5Hz-V','5.5Hz-AV','5.5Hz-A','40Hz-V','40Hz-AV','40Hz-A','80Hz-V','80Hz-AV','80Hz-A'}; %define conditions
flickerfreq_conditions={'5.5Hz','8Hz','11Hz','14Hz','17Hz','20Hz','23Hz','26Hz','29Hz','32Hz','35Hz','38Hz','40Hz','42Hz','45Hz','48Hz','51Hz','54Hz','57Hz','63Hz','66Hz','69Hz','72Hz','75Hz','78Hz','80Hz'};

sessions=struct;
examples=cell([0 6]);
for s=1:size(all_sessions,1)
    disp(['Loading data for session ' num2str(s) ' out of ' num2str(size(all_sessions,1))]);
    sessions(s).subjectID=all_sessions{s,'sub'};
    sessions(s).task=all_sessions{s,'task'};
    sessions(s).ses=all_sessions{s,'ses'};
    sessions(s).modality=all_sessions{s,'modality'};

    %load modulation data:
    sessions(s).ssep_sig=readtable([root_dir '/stg-analyses/task-' sessions(s).task{:} '/sub-' sessions(s).subjectID{:} '/ses-' sessions(s).ses{:} '/LFP/static_ent/LFP_pvalue_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
    sessions(s).ssep_amp=readtable([root_dir '/stg-analyses/task-' sessions(s).task{:} '/sub-' sessions(s).subjectID{:} '/ses-' sessions(s).ses{:} '/LFP/static_ent/LFP_zscore_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
    sessions(s).anat=importdata([root_dir '/stg-preproc/sub-' sessions(s).subjectID{:} '/anat/electrode-models/sub-' sessions(s).subjectID{:} '_electrodes_info.mat']);
    
    sessions(s).data=importdata([root_dir '/stg-analyses/task-' sessions(s).task{:} '/sub-' sessions(s).subjectID{:} '/ses-' sessions(s).ses{:} '/LFP/static_ent/'...
                                    'sub-' sessions(s).subjectID{:} '_stg-analysis_task-' sessions(s).task{:} '_ses-' sessions(s).ses{:} '_nat-echo-analysis-preproc_cleaned-linenoise.mat'],'echo_results');

    %select subset of scenarios where: 1) have signiticant modulation and 2) ABC shows
    %significant cycles for at x percent of cycles:
    sessions(s).data_subset=sessions(s).data;
    if strcmp(sessions(s).task,'flickerfreq')
        temp_conditions=strcat(flickerfreq_conditions,'-',sessions(s).modality);
    else
        temp_conditions=flickerneuro_conditions;
    end
    for c=1:length(temp_conditions) %for each condition
        chs_interest=zeros([1 length(sessions(s).data{1}.label)]);
        bin_times=(-ceil(str2double(regexprep(temp_conditions{c},'Hz-.+','')))*1/(str2double(regexprep(temp_conditions{c},'Hz-.+',''))):1/str2double(regexprep(temp_conditions{c},'Hz-.+','')):11);
        bin_times(end)=[];
        for ch=1:length(sessions(s).data{1}.label) %for each channel
            if sessions(s).ssep_sig{sessions(s).data{1}.label{ch},temp_conditions{c}}<0.05 %if we have significant modulation
                sessions(s).percent_cycles(ch,c)=sum(sessions(s).data{c}.evoked_hil_rec_bin(bin_times>=0,ch)>sessions(s).data{c}.threshold)/(str2double(regexprep(temp_conditions{c},'Hz-.+',''))*10);
                temp1=sessions(s).data{c}.evoked_hil_rec_bin(bin_times>=10,ch)>sessions(s).data{c}.threshold;
                temp2=0;
                for i=1:length(temp1)
                    if temp1(i)
                        temp2=temp2+1;
                    else
                        break;
                    end
                end
                sessions(s).echo_pers(ch,c)=temp2;
                if sessions(s).percent_cycles(ch,c)>0.1 %if we have at least 1/10th of cycles during stim that are significant
                    %calculate number of persistent cycles and save this
                    %channel as of interest:
                    chs_interest(ch)=1;
                    examples(end+1,:)={s temp_conditions{c} sessions(s).data{1}.label{ch} temp2 sessions(s).ssep_amp{sessions(s).data{1}.label{ch},temp_conditions{c}} sum(sessions(s).data{c}.evoked_hil_rec_bin(bin_times>=0 & bin_times<10,ch)>sessions(s).data{c}.threshold)/(str2double(regexprep(temp_conditions{c},'Hz-.+',''))*10)};
                end
            elseif sessions(s).ssep_sig{sessions(s).data{1}.label{ch},temp_conditions{c}}>=0.05 %if we don't have significant modulation
                sessions(s).percent_cycles(ch,c)=-1;
                sessions(s).echo_pers(ch,c)=-1;
            end
        end
        
        sessions(s).data_subset{c}.label(~chs_interest)=[];
        sessions(s).data_subset{c}.evoked(~chs_interest,:)=[];
        sessions(s).data_subset{c}.evoked_filtered(~chs_interest,:)=[];
        sessions(s).data_subset{c}.evoked_hil(~chs_interest,:)=[];
        sessions(s).data_subset{c}.evoked_hil_z(~chs_interest,:)=[];
        sessions(s).data_subset{c}.evoked_hil_rec(~chs_interest,:)=[];
        sessions(s).data_subset{c}.evoked_hil_rec_bin(:,~chs_interest)=[];
        sessions(s).data_subset{c}.good_channels(~chs_interest)=[];
        sessions(s).data_subset{c}.n_cycles(~chs_interest)=[];
        sessions(s).data_subset{c}.onsets(~chs_interest)=[];
    end
    sessions(s).percent_cycles=array2table(sessions(s).percent_cycles,'VariableNames',temp_conditions,'RowNames',sessions(s).data{1}.label);
    sessions(s).echo_pers=array2table(sessions(s).echo_pers,'VariableNames',temp_conditions,'RowNames',sessions(s).data{1}.label);
     
 
    for c=1:length(sessions(s).data)
        sessions(s).data{c}=rmfield(sessions(s).data{c},{'evoked','evoked_filtered','evoked_hil','evoked_hil_rec'});
    end
end


%% plot for figure:

%plot 1 example full scale:
time_window_start=[-0.1 0.5]; %time window of start of trial that want zoom-in
time_window_end=[9.8 10.3]; %time window of end of trial that want zoom-in
for i=1
    s=find(strcmp([sessions.subjectID],selected_examples{i,1}) & strcmp([sessions.task],selected_examples{i,2}) & strcmp([sessions.ses],selected_examples{i,3}));
    if strcmp(sessions(s).task,'flickerfreq')
        temp_conditions=strcat(flickerfreq_conditions,'-',sessions(s).modality);
    else
        temp_conditions=flickerneuro_conditions;
    end
    c=selected_examples{i,4};
    ch=find(strcmp(sessions(s).data_subset{strcmp(temp_conditions,c)}.label,selected_examples{i,5}));
    
    %define times, at sampling rate:
    times=-1:1/sessions(s).data{strcmp(temp_conditions,c)}.sampleRate:11;
    times(end)=[];

    %define mid-bin times:
    bin_times=(-ceil(str2double(regexprep(c,'Hz-.+','')))*1/(str2double(regexprep(c,'Hz-.+',''))):1/str2double(regexprep(c,'Hz-.+','')):11)+1/str2double(regexprep(c,'Hz-.+',''))/2;
    bin_times(end)=[];

    %define indices of times at zoom-in start and end:
    t_index_start=times>=time_window_start(1) & times<=time_window_start(2);
    t_index_end=times>=time_window_end(1) & times<=time_window_end(2);

    %define bins that were significant:
    bins_of_interest=sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked_hil_rec_bin(:,ch)>sessions(s).data_subset{strcmp(temp_conditions,c)}.threshold;

    %define bins that we are going to highlight in zoomed-in plots:
    temp=bin_times(bins_of_interest);
    bin_times_start=temp(temp>=time_window_start(1) & temp<=time_window_start(2));
    bin_times_end=temp(temp>=time_window_end(1) & temp<=time_window_end(2));
    
    %figure('Name',['Subject ' sessions(s).subjectID{:} ', condition ' c ', channel index ' num2str(ch)],'Units','normalized','Position',[0.1 0.1 0.8 0.8]);
    figure_making('width',5.2,'height',2.5);
    tiledlayout(5,2,'Padding','none','TileSpacing','none');

    %plot evoked potential:
    nexttile(1,[1 2]);
    plot(times,sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked(ch,:),'k');
    xline(0,'Color','k');
    xline(10,'Color','k');
    xline(sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),'--','Color','k');
    xline(10+sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),'--','Color','k');
    %title('Evoked response (averaged across trials)');
    temp=gca;
    %text(-1.6,temp.YLim(2)+range(temp.YLim)*0.15,[sessions(s).subjectID{:} ', condition ' c ', channel ' selected_examples{i,3}],'FontWeight','bold');
    %xlabel('Time (s)');
    plot_flicker(str2double(regexprep(c,'Hz-.+','')),temp.YLim,10,condition_color(c),0);
    axis tight;
    xlim([-1 11]);
    ylabel('Potential (\muV)')
    box off;
    ax=gca;
    ax.YAxis.Visible='off';
    ax.XTick=[0 10];
    xline(time_window_start(1),'Color',[0.5 0.5 0.5]);
    xline(time_window_start(2),'Color',[0.5 0.5 0.5]);
    xline(time_window_end(1),'Color',[0.5 0.5 0.5]);
    xline(time_window_end(2),'Color',[0.5 0.5 0.5]);

    %plot zoom-in of first and last 500ms of evoked potential:
    y_mean_start=mean(sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked(ch,t_index_start));
    y_half_range_start=max(abs(y_mean_start-sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked(ch,t_index_start)));
    y_mean_end=mean(sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked(ch,t_index_end));
    y_half_range_end=max(abs(y_mean_end-sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked(ch,t_index_end)));
    
    if y_half_range_start>y_half_range_end
        y_half_range=y_half_range_start;
    else
        y_half_range=y_half_range_end;
    end
    
    nexttile(3);
    plot(times(t_index_start),sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked(ch,t_index_start),'k');
    hold on;
    xline(0,'Color','k');
    xline(sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),'--','Color','k');
    scatter(bin_times_start+sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),repmat(y_mean_start+y_half_range+abs(y_mean_start-y_half_range)*0.3,1,length(bin_times_start)),'r','filled','SizeData',5);
    xlim([time_window_start(1) time_window_start(2)]);
    ylim([y_mean_start-y_half_range y_mean_start+y_half_range+abs(y_mean_start-y_half_range)*0.6]);
    temp=gca;
    plot_flicker(str2double(regexprep(c,'Hz-.+','')),temp.YLim,time_window_start(2),condition_color(c),0);
    ylim([y_mean_start-y_half_range y_mean_start+y_half_range+abs(y_mean_start-y_half_range)*0.6+y_half_range/3]);
    box off;
    ax=gca;
    ax.YAxis.Visible='off';
    ax.XAxis.Visible='off';
    
    nexttile(4);
    plot(times(t_index_end),sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked(ch,t_index_end),'k');
    hold on;
    xline(10,'Color','k');
    xline(10+sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),'--','Color','k');
    scatter(bin_times_end+sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),repmat(y_mean_end+y_half_range+abs(y_mean_end-y_half_range)*0.3,1,length(bin_times_end)),'r','filled','SizeData',5);
    xlim([time_window_end(1) time_window_end(2)]);
    ylim([y_mean_end-y_half_range y_mean_end+y_half_range+abs(y_mean_start-y_half_range)*0.6]);
    temp=gca;
    plot_flicker(str2double(regexprep(c,'Hz-.+','')),temp.YLim,10,condition_color(c),0);
    ylim([y_mean_start-y_half_range y_mean_start+y_half_range+abs(y_mean_start-y_half_range)*0.6+y_half_range/3]);
    box off;
    ax=gca;
    ax.YAxis.Visible='off';
    ax.XAxis.Visible='off';
    

    %plot zoom-in of first and last 500ms of filtered potential:
    y_mean_start=mean(sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked_filtered(ch,t_index_start));
    y_half_range_start=max(abs(y_mean_start-sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked_filtered(ch,t_index_start)));
    y_mean_end=mean(sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked_filtered(ch,t_index_end));
    y_half_range_end=max(abs(y_mean_start-sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked_filtered(ch,t_index_end)));
    
    if y_half_range_start>y_half_range_end
        y_half_range=y_half_range_start;
    else
        y_half_range=y_half_range_end;
    end

    nexttile(5);
    plot(times(t_index_start),sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked_filtered(ch,t_index_start),'k');
    hold on;
    xline(0,'Color','k');
    xline(sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),'--','Color','k');
    temp_t_index=[];
    for t=bin_times_start+sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch)
        [~,temp_min]=min(abs(times-t));
        temp_t_index=[temp_t_index temp_min];
    end
    scatter(bin_times_start+sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),repmat(y_mean_start+y_half_range+abs(y_mean_start-y_half_range)*0.3,1,length(bin_times_start)),'r','filled','SizeData',5);
    xlim([time_window_start(1) time_window_start(2)]);
    ylim([y_mean_start-y_half_range y_mean_start+y_half_range+abs(y_mean_start-y_half_range)*0.6]);
    temp=gca;
    plot_flicker(str2double(regexprep(c,'Hz-.+','')),temp.YLim,time_window_start(2),condition_color(c),0);
    ylim([y_mean_start-y_half_range y_mean_start+y_half_range+abs(y_mean_start-y_half_range)*0.6+y_half_range/3]);
    box off;
    ax=gca;
    ax.YAxis.Visible='off';
    ax.XTick=[0 0.4];
    
    nexttile(6);
    plot(times(t_index_end),sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked_filtered(ch,t_index_end),'k');
    hold on;
    xline(10,'Color','k');
    xline(10+sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),'--','Color','k');
    scatter(bin_times_end+sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),repmat(y_mean_end+y_half_range+abs(y_mean_end-y_half_range)*0.3,1,length(bin_times_end)),'r','filled','SizeData',5);
    xlim([time_window_end(1) time_window_end(2)]);
    ylim([y_mean_end-y_half_range y_mean_end+y_half_range+abs(y_mean_start-y_half_range)*0.6]);
    temp=gca;
    plot_flicker(str2double(regexprep(c,'Hz-.+','')),temp.YLim,10,condition_color(c),0);
    ylim([y_mean_start-y_half_range y_mean_start+y_half_range+abs(y_mean_start-y_half_range)*0.6+y_half_range/3]);
    box off;
    ax=gca;
    ax.YAxis.Visible='off';
    ax.XTick=[9.85 10];
    
    
    %plot filtered signal:
    nexttile(7,[1 2]);
    plot(times,sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked_filtered(ch,:),'k');
    %yline(0,'--');
    xline(0,'Color','k');
    xline(10,'Color','k');
    xline(sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),'--','Color','k');
    xline(10+sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),'--','Color','k');
    xlim([-1 11]);
    temp=gca;
    plot_flicker(str2double(regexprep(c,'Hz-.+','')),temp.YLim,10,condition_color(c),0);
    box off;
    ax=gca;
    ax.YAxis.Visible='off';
    ax.XTick=[0 10];
    xline(time_window_start(1),'Color',[0.5 0.5 0.5]);
    xline(time_window_start(2),'Color',[0.5 0.5 0.5]);
    xline(time_window_end(1),'Color',[0.5 0.5 0.5]);
    xline(time_window_end(2),'Color',[0.5 0.5 0.5]);
    
    %plot time-frequency plot to show specificity of lingering response at frequency of stimulation:
    nexttile(9,[1 2]);
    low_freq_lim=str2double(regexprep(c,'Hz-.+',''))-30;
    if low_freq_lim<0
        low_freq_lim=0;
    end
    [temp_wt,f]=cwt(sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked(ch,:),sessions(s).data{strcmp(temp_conditions,c)}.sampleRate,'FrequencyLimits',[low_freq_lim str2double(regexprep(c,'Hz-.+',''))+30]);
    pcolor(times,f,abs(temp_wt));
    shading interp;
    xline(0,'Color','k');
    xline(10,'Color','k');
    xline(sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),'--','Color','k');
    xline(10+sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),'--','Color','k');
    %yline(str2double(regexprep(c,'Hz-.+','')),'--','Color',[0.5 0.5 0.5]);
    ax=gca;
    ax.YTick=str2num(regexprep(c,'Hz-.+',''));
    ax.YTick=[low_freq_lim+3 ax.YTick ax.YTick+27];
    %ax.ColorScale='log';
    xlabel('Time (s)');

    figure_making('operation','save','filename',[figures_dir '/Figure_S8A.svg']);

end


%plot 6 examples zoomed-in:
figure_making('width',6.5,'height',3.5);
tiledlayout(6,4,'Padding','none','TileSpacing','none');
time_window_start=[-0.1 0.5]; %time window of start of trial that want zoom-in
time_window_end=[9.8 10.3]; %time window of end of trial that want zoom-in
% time_window_start=[-0.05 0.15]; %time window of start of trial that want zoom-in
% time_window_end=[9.95 10.15]; %time window of end of trial that want zoom-in
subplot_nber=1;
for i=[1 4 2 5 3 6]
    s=find(strcmp([sessions.subjectID],selected_examples{i,1}) & strcmp([sessions.task],selected_examples{i,2}) & strcmp([sessions.ses],selected_examples{i,3}));
    if strcmp(sessions(s).task,'flickerfreq')
        temp_conditions=strcat(flickerfreq_conditions,'-',sessions(s).modality);
    else
        temp_conditions=flickerneuro_conditions;
    end
    c=selected_examples{i,4};
    ch=find(strcmp(sessions(s).data_subset{strcmp(temp_conditions,c)}.label,selected_examples{i,5}));
    
    %define times, at sampling rate:
    times=-1:1/sessions(s).data{strcmp(temp_conditions,c)}.sampleRate:11;
    times(end)=[];

    %define mid-bin times:
    bin_times=(-ceil(str2double(regexprep(c,'Hz-.+','')))*1/(str2double(regexprep(c,'Hz-.+',''))):1/str2double(regexprep(c,'Hz-.+','')):11)+1/str2double(regexprep(c,'Hz-.+',''))/2;
    bin_times(end)=[];

    %define indices of times at zoom-in start and end:
    t_index_start=times>=time_window_start(1) & times<=time_window_start(2);
    t_index_end=times>=time_window_end(1) & times<=time_window_end(2);

    %define bins that were significant:
    bins_of_interest=sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked_hil_rec_bin(:,ch)>sessions(s).data_subset{strcmp(temp_conditions,c)}.threshold;

    %define bins that we are going to highlight in zoomed-in plots:
    temp=bin_times(bins_of_interest);
    bin_times_start=temp(temp>=time_window_start(1) & temp<=time_window_start(2));
    bin_times_end=temp(temp>=time_window_end(1) & temp<=time_window_end(2));
    
    %figure('Name',['Subject ' sessions(s).subjectID{:} ', condition ' c ', channel index ' num2str(ch)],'Units','normalized','Position',[0.1 0.1 0.8 0.8]);
    

    %plot zoom-in of first and last 500ms of evoked potential:
    y_mean_start=mean(sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked(ch,t_index_start));
    y_half_range_start=max(abs(y_mean_start-sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked(ch,t_index_start)));
    y_mean_end=mean(sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked(ch,t_index_end));
    y_half_range_end=max(abs(y_mean_end-sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked(ch,t_index_end)));
    
    if y_half_range_start>y_half_range_end
        y_half_range=y_half_range_start;
    else
        y_half_range=y_half_range_end;
    end
    
    nexttile(subplot_nber);
    plot(times(t_index_start),sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked(ch,t_index_start),'k');
    hold on;
    xline(0,'Color','k');
    xline(sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),'--','Color','k');
    scatter(bin_times_start+sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),repmat(y_mean_start+y_half_range+abs(y_mean_start-y_half_range)*0.3,1,length(bin_times_start)),'r','filled','SizeData',4);
    ylim([y_mean_start-y_half_range y_mean_start+y_half_range+abs(y_mean_start-y_half_range)*0.6]);
    temp=gca;
    plot_flicker(str2double(regexprep(c,'Hz-.+','')),temp.YLim,time_window_start(2)+1,condition_color(c),0);
    ylim([y_mean_start-y_half_range y_mean_start+y_half_range+abs(y_mean_start-y_half_range)*0.6+y_half_range/3]);
    xlim([time_window_start(1) time_window_start(2)]);
    box off;
    ax=gca;
    ax.YAxis.Visible='off';
    ax.XAxis.Visible='off';
    
    nexttile(subplot_nber+1);
    plot(times(t_index_end),sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked(ch,t_index_end),'k');
    hold on;
    xline(10,'Color','k');
    xline(10+sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),'--','Color','k');
    scatter(bin_times_end+sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),repmat(y_mean_end+y_half_range+abs(y_mean_end-y_half_range)*0.3,1,length(bin_times_end)),'r','filled','SizeData',4);
    ylim([y_mean_end-y_half_range y_mean_end+y_half_range+abs(y_mean_start-y_half_range)*0.6]);
    temp=gca;
    plot_flicker(str2double(regexprep(c,'Hz-.+','')),temp.YLim,10,condition_color(c),0);
    ylim([y_mean_end-y_half_range y_mean_end+y_half_range+abs(y_mean_start-y_half_range)*0.6+y_half_range/3]);
    xlim([time_window_end(1) time_window_end(2)]);
    box off;
    ax=gca;
    ax.YAxis.Visible='off';
    ax.XAxis.Visible='off';
    
    nexttile(subplot_nber+4);
    low_freq_lim=str2double(regexprep(c,'Hz-.+',''))-30;
    if low_freq_lim<0
        low_freq_lim=0;
    end
    [temp_wt,f]=cwt(sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked(ch,t_index_start),sessions(s).data{strcmp(temp_conditions,c)}.sampleRate,'FrequencyLimits',[low_freq_lim str2double(regexprep(c,'Hz-.+',''))+30]);
    pcolor(times(t_index_start),f,abs(temp_wt));
    shading interp;
    xline(0,'Color','k');
    xline(sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),'--','Color','k');
    ax=gca;
    ax.YTick=str2num(regexprep(c,'Hz-.+',''));
    ax.YTick=[low_freq_lim+3 ax.YTick ax.YTick+27];
    if subplot_nber+4~=21
        ax.XAxis.Visible='off';
    end
    
    
    nexttile(subplot_nber+5);
    [temp_wt,f]=cwt(sessions(s).data_subset{strcmp(temp_conditions,c)}.evoked(ch,t_index_end),sessions(s).data{strcmp(temp_conditions,c)}.sampleRate,'FrequencyLimits',[low_freq_lim str2double(regexprep(c,'Hz-.+',''))+30]);
    pcolor(times(t_index_end),f,abs(temp_wt));
    shading interp;
    xline(0,'Color','k');
    xline(sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),'--','Color','k');
    xline(10,'Color','k');
    xline(10+sessions(s).data_subset{strcmp(temp_conditions,c)}.onsets(ch),'--','Color','k');
    ax=gca;
    if subplot_nber+5~=22
        ax.XAxis.Visible='off';
    end
    ax.YAxis.Visible='off';
    
    subplot_nber=subplot_nber+2;
    
    if mod(subplot_nber,4)==1
        subplot_nber=subplot_nber+4;
    end
end

figure_making('operation','save','filename',[figures_dir '/Figure_S8B.svg']);

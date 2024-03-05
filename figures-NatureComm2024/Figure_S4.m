%Parts of figure: S4.
%2024/02/26

%% set directories:

root_dir=define_flicker_root_dir;
source_data=[root_dir '/stg-analyses/NatureComm2024-figures/Source Data']; %where source data for figures is stored
figures_dir=[root_dir '/stg-analyses/NatureComm2024-figures'];

%% fetch data:
[~,all_sessions]=fetch_flicker_subjectIDs([define_flicker_root_dir],'flickerneuro');

sessions=struct();
for s=1:size(all_sessions,1)
    
    disp(['Loading and preprocessing session ' num2str(s) ' out of ' num2str(size(all_sessions,1)) '...']);
    
    %assign some information:
    sessions(s).subjectID=all_sessions{s,'sub'}{:};
    sessions(s).task=all_sessions{s,'task'}{:};
    sessions(s).ses=all_sessions{s,'ses'}{:};
    
    %fetch significance and amplitude of modulation for all channels and
    %conditions for this session:
    temp_sigtable=readtable([root_dir '/stg-analyses/task-' all_sessions{s,'task'}{:} '/sub-' all_sessions{s,'sub'}{:} '/ses-' all_sessions{s,'ses'}{:} '/LFP/static_ent/LFP_pvalue_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
    temp_amptable=readtable([root_dir '/stg-analyses/task-' all_sessions{s,'task'}{:} '/sub-' all_sessions{s,'sub'}{:} '/ses-' all_sessions{s,'ses'}{:} '/LFP/static_ent/LFP_zscore_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);

    %initialize tables that will keep track of significance, amplitude or
    %all channels for periodic and rand conditions, response at 30Hz, 40Hz
    %and 50Hz.
    sessions(s).sig_40Hz=array2table(nan([size(temp_sigtable,1),9]),'VariableNames',{'V-30Hz-resp','V-40Hz-resp','V-50Hz-resp','AV-30Hz-resp','AV-40Hz-resp','AV-50Hz-resp','A-30Hz-resp','A-40Hz-resp','A-50Hz-resp'},'RowNames',temp_sigtable.Properties.RowNames);
    sessions(s).amp_40Hz=array2table(nan([size(temp_sigtable,1),9]),'VariableNames',{'V-30Hz-resp','V-40Hz-resp','V-50Hz-resp','AV-30Hz-resp','AV-40Hz-resp','AV-50Hz-resp','A-30Hz-resp','A-40Hz-resp','A-50Hz-resp'},'RowNames',temp_sigtable.Properties.RowNames);

    sessions(s).sig_R=array2table(nan([size(temp_sigtable,1),9]),'VariableNames',{'V-30Hz-resp','V-40Hz-resp','V-50Hz-resp','AV-30Hz-resp','AV-40Hz-resp','AV-50Hz-resp','A-30Hz-resp','A-40Hz-resp','A-50Hz-resp'},'RowNames',temp_sigtable.Properties.RowNames);
    sessions(s).amp_R=array2table(nan([size(temp_sigtable,1),9]),'VariableNames',{'V-30Hz-resp','V-40Hz-resp','V-50Hz-resp','AV-30Hz-resp','AV-40Hz-resp','AV-50Hz-resp','A-30Hz-resp','A-40Hz-resp','A-50Hz-resp'},'RowNames',temp_sigtable.Properties.RowNames);

    sessions(s).sig_40Hz(:,{'V-40Hz-resp','AV-40Hz-resp','A-40Hz-resp'})=temp_sigtable(:,{'40Hz-V','40Hz-AV','40Hz-A'});
    sessions(s).amp_40Hz(:,{'V-40Hz-resp','AV-40Hz-resp','A-40Hz-resp'})=temp_amptable(:,{'40Hz-V','40Hz-AV','40Hz-A'});

    %load PSDs:
    temp_psd=importdata([root_dir '/stg-preproc/sub-' all_sessions{s,'sub'}{:} '/task-' all_sessions{s,'task'}{:} '/ses-' all_sessions{s,'ses'}{:} '/LFP/static_ent/sub-' all_sessions{s,'sub'}{:} '_stg-analysis_task-' all_sessions{s,'task'}{:} '_ses-' all_sessions{s,'ses'}{:} '_nat-psd-refLaplacian.mat']);

    %complete tables:
    control_condition='Baseline';
    
    %THINK ABOUT DOING CALCULATIONS ONLY FOR THOSE MODULATIONS THAT ARE
    %SIGNIFICANT AT 40HZ?
    for nat={'40Hz','R'} %for periodic conditions, and random conditions
        for mod={'V','AV','A'} %for each modality
            for freq=[30 40 50] %for response at 30Hz, 40Hz and 50Hz
                for ch=1:size(sessions(s).sig_40Hz,1) %for each channel

                    %check that channel labels match:
                    if ~strcmp(sessions(s).sig_40Hz.Properties.RowNames{ch},temp_psd.label{ch})
                        error('Channel labels do not match');
                    end
                    
                    if ~isequal([nat freq],{'40Hz',40}) %already have values if this is periodic condition, response at 40Hz
                        %initialize stim and baseline values:
                        stim_values=[];
                        baseline_values=[];
                        
                        [~,freq_interest_index]=min(abs(temp_psd.data{1,strcmp(temp_psd.condition,[nat{:} '-' mod{:}])}{3}-freq)); %find index of frequency closest to frequency of interest- have to do this because sample rate of EDF file not an integer sometimes (error with Natus)
                        
                        for tr=1:size(temp_psd.data{ch,strcmp(temp_psd.condition,[nat{:} '-' mod{:}])}{1},1) %for however many number of trials of given condition there are
                            current_stim_value=temp_psd.data{ch,strcmp(temp_psd.condition,[nat{:} '-' mod{:}])}{1}(tr,freq_interest_index);
                            current_baseline_value=temp_psd.data{ch,strcmp(temp_psd.condition,control_condition)}{1}(tr,freq_interest_index);
            
                            stim_values=[stim_values current_stim_value];
                            baseline_values=[baseline_values current_baseline_value];
                        end
    
                        sessions(s).(['sig_' nat{:}]){ch,[mod{:} '-' num2str(freq) 'Hz-resp']}=pval_randomshuffle([stim_values' baseline_values'],10000);
                        sessions(s).(['amp_' nat{:}]){ch,[mod{:} '-' num2str(freq) 'Hz-resp']}=(mean(stim_values)/mean(baseline_values))-1;
                    end
                end
            end
        end
    end

end

save([root_dir '/stg-analyses/task-flickerneuro/mod40Hz_40HzvsRstim.mat'],'sessions','-v7.3');

%% aggregate and summarize results:

sessions=importdata([root_dir '/stg-analyses/task-flickerneuro/mod40Hz_40HzvsRstim.mat']);

%take the difference between modulation at 40Hz vs 30Hz and 40Hz vs 50Hz:
mat_40Hz_3050=[];
mat_R_3050=[];
scenarios_of_interest={};
sessions_included=[];
for s=1:length(sessions)
    for ch=1:size(sessions(s).sig_40Hz,1)
        for mod={'V','AV','A'}
            if sessions(s).sig_40Hz{ch,[mod{:} '-40Hz-resp']}<0.05 %if we have significant entrainment for that channel and condition at 40Hz
                mat_40Hz_3050(end+1)=sessions(s).amp_40Hz{ch,[mod{:} '-40Hz-resp']}-mean([sessions(s).amp_40Hz{ch,[mod{:} '-30Hz-resp']} sessions(s).amp_40Hz{ch,[mod{:} '-50Hz-resp']}]);
                scenarios_of_interest(end+1,:)={sessions(s).subjectID sessions(s).sig_40Hz.Properties.RowNames{ch} mod{:} sessions(s).amp_40Hz{ch,[mod{:} '-40Hz-resp']} mat_40Hz_3050(end)};
                sessions_included=[sessions_included s];
            end
            

            if sessions(s).sig_R{ch,[mod{:} '-40Hz-resp']}<0.05 %if we have significant entrainment for that channel and condition at 40Hz
                mat_R_3050(end+1)=sessions(s).amp_R{ch,[mod{:} '-40Hz-resp']}-mean([sessions(s).amp_R{ch,[mod{:} '-30Hz-resp']} sessions(s).amp_R{ch,[mod{:} '-50Hz-resp']}]);
                %scenarios_of_interest(end+1,:)={sessions(s).subjectID sessions(s).sig_40Hz.Properties.RowNames{ch} mod{:} sessions(s).amp_R{ch,[mod{:} '-40Hz-resp']} mat_R_3050(end)};
                sessions_included=[sessions_included s];
            end
        end
    end
end

disp(['Total number of sessions included: ' num2str(length(unique(sessions_included)))]);

%plot distributions of differences:
hist_tbl=table('Size',[21 3],'VariableTypes',{'string','double','double'},'VariableNames',{'bins','percent_ch_mod_40Hz_stim','percent_ch_mod_40Hz_rand'});
hist_edges=[-Inf 0:0.1:1 2:10];
hist_tbl(:,1)={'<=0';'>0-0.1';'0.1-0.2';'0.2-0.3';'0.3-0.4';'0.4-0.5';'0.5-0.6';'0.6-0.7';'0.7-0.8';'0.8-0.9';'0.9-1';'1-2';'2-3';'3-4';'4-5';'5-6';'6-7';'7-8';'8-9';'9-10';'>10'};
hist_mat=nan(1,length(hist_edges));
for e=1:length(hist_edges)-1
    hist_mat(e)=sum(mat_40Hz_3050>hist_edges(e) & mat_40Hz_3050<=hist_edges(e+1));
end
hist_mat(end)=sum(mat_40Hz_3050>hist_edges(end));
hist_mat=hist_mat/sum(hist_mat)*100;
hist_tbl{:,2}=hist_mat';
hist_mat=[hist_mat(1) 0 hist_mat(2:11) 0 hist_mat(12:20) 0 hist_mat(end)];
figure_making('width',6.5,'height',2);
histogram('BinEdges',0:length(hist_edges)+3,'BinCounts',hist_mat,'FaceColor','r');
hold on;
hist_mat=zeros(1,length(hist_edges));
for e=1:length(hist_edges)-1
    hist_mat(e)=sum(mat_R_3050>hist_edges(e) & mat_R_3050<=hist_edges(e+1));
end
hist_mat(end)=sum(mat_R_3050>hist_edges(end));
hist_mat=hist_mat/sum(hist_mat)*100;
hist_tbl{:,3}=hist_mat';
hist_mat=[hist_mat(1) 0 hist_mat(2:11) 0 hist_mat(12:20) 0 hist_mat(end)];
histogram('BinEdges',0:length(hist_edges)+3,'BinCounts',hist_mat,'FaceColor',[0.5 0.5 0.5]);
ax=gca;
ax.XTick=[0.5 2:12 13:22 23.5];
ax.XTickLabel=[{'\leq0'} {'>0'} regexprep(cellstr(num2str([0.1:0.1:1]')),' ','')' regexprep(cellstr(num2str([1:10]')),' ','')' {'>10'}];
xline(1.5,'--','Color','k');
xline(12.5,'--','Color','k');
xline(22.5,'--','Color','k');
ylabel('% of channels modulated at 40Hz');
xlabel('Modulation amplitude difference between 40Hz and 30/50Hz');
box off;
text(23,ax.YLim(end)-range(ax.YLim)/12,['\color[rgb]{' num2str([1 0 0]) '}40Hz stim'],'FontSize',7);
text(23,ax.YLim(end)-(range(ax.YLim)/12)*2,['\color[rgb]{' num2str([0.5 0.5 0.5]) '}R stim'],'FontSize',7);

%show number of significant contacts for 40Hz stim and R stim:


%save figure and associated source data:
figure_making('operation','save','filename',[figures_dir '/Figure_S4C.pdf']);
writetable(hist_tbl,[source_data '/Figure_S4C.csv']);


%% plot examples:

%examples modulation in sensory regions:
examples=readtable([figures_dir '/NatureComm2024_examples.xlsx']);
examples=table2array(examples(startsWith(examples.example,'Figure_S4'),:));
[~,temp_index]=sort(examples(:,1));
examples=examples(temp_index,[2 5 6 1]);

%fetch data for those subjects:
for subjectID=unique(examples(:,1))'
    ses_nber=find(strcmp({sessions.subjectID},subjectID));
    sessions(ses_nber).psd_data=importdata([define_flicker_root_dir '/stg-preproc/sub-' sessions(ses_nber).subjectID '/task-' sessions(ses_nber).task '/ses-' sessions(ses_nber).ses '/LFP/static_ent/sub-' sessions(ses_nber).subjectID '_stg-analysis_task-' sessions(ses_nber).task '_ses-' sessions(ses_nber).ses '_nat-psd-refLaplacian.mat']);
    if strcmp(subjectID,examples{1,1})
        sessions(ses_nber).lfp_data=importdata([define_flicker_root_dir '/stg-preproc/sub-' sessions(ses_nber).subjectID '/task-' sessions(ses_nber).task '/ses-' sessions(ses_nber).ses '/LFP/static_ent/sub-' sessions(ses_nber).subjectID '_stg-preproc_task-' sessions(ses_nber).task '_ses-' sessions(ses_nber).ses '_nat-trialsdata-refLaplacian-preproc.mat']);
        sessions(ses_nber).trials_data=importdata([define_flicker_root_dir '/stg-preproc/sub-' sessions(ses_nber).subjectID '/task-' sessions(ses_nber).task '/ses-' sessions(ses_nber).ses '/sub-' sessions(ses_nber).subjectID '_stg-preproc_task-' sessions(ses_nber).task '_ses-' sessions(ses_nber).ses '_nat-beh.mat']);
    end
end

[~,repo]=define_flicker_root_dir;
addpath(genpath(repo.fieldtrip)); %code below uses the FieldTrip's nearest function, as opposed to the MATLAB's; hence, puting the FieldTrip repo to top of path

%plot PSDs:
for ex=1:3
    mod_interest=regexprep(examples{ex,3},'.+-','');
    freq_interest=regexprep(examples{ex,3},'Hz-.+','');

    if ex==1 %plot example trial for 40Hz and random conditions
        figure_making('width',6.5,'height',1);
        tiledlayout(1,2,'Padding','compact');
        ses_nber=find(strcmp({sessions.subjectID},examples{ex,1}));
        
        
        nexttile(1);
        trials_of_interest=find(sessions(ses_nber).lfp_data.trialinfo==find(strcmp(sessions(ses_nber).trials_data.condition_code,examples{ex,3})));
        cfg=[];
        cfg.trials=trials_of_interest;
        cfg.keeptrials='no';
        temp=ft_timelockanalysis(cfg,sessions(ses_nber).lfp_data);

        plot(temp.time,temp.avg(strcmp(temp.label,examples{ex,2}),:),'k');
        xlim([-0.1 0.4]);
        temp_ylim=gca;
        temp_ylim=[min(temp_ylim.Children.YData(temp_ylim.Children.XData>=-0.1 & temp_ylim.Children.XData<=0.5)) max(temp_ylim.Children.YData(temp_ylim.Children.XData>=-0.1 & temp_ylim.Children.XData<=0.5))];
        first_plot_ylim=get(gca,'ylim');
        plot_flicker(str2double(freq_interest),temp_ylim,10,condition_color(examples{ex,3}),0);
        first_plot_ylim2=ylim;
        %ylim([temp(1) temp(2)+((temp(2)-temp(1))/5)]);
        box off;
        xline(0,'--','LineWidth',1,'Label','','LabelHorizontalAlignment','left');
        set(gca,'visible','off');
        ax=gca;
        annotation('line',[ax.Position(1) ax.Position(1)+0.1*(ax.Position(3)/range(ax.XLim))],[ax.Position(2) ax.Position(2)],'Linewidth',1) %bar represents 0.1s
        annotation('line',[ax.Position(1) ax.Position(1)],[ax.Position(2) ax.Position(2)+50*(ax.Position(4)/range(ax.YLim))],'Linewidth',1) %bar represents 50microV
        
        nexttile(2);
        trials_of_interest=find(sessions(ses_nber).lfp_data.trialinfo==find(strcmp(sessions(ses_nber).trials_data.condition_code,['R-' mod_interest])));
        cfg=[];
        cfg.trials=trials_of_interest;
        cfg.keeptrials='no';
        temp=ft_timelockanalysis(cfg,sessions(ses_nber).lfp_data);

        plot(temp.time,temp.avg(strcmp(temp.label,examples{ex,2}),:),'k');
        %plot(sessions(ses_nber).lfp_data.time{1},sessions(ses_nber).lfp_data.trial{trials_of_interest(1)}(strcmp(sessions(ses_nber).lfp_data.label,examples{ex,2}),:),'k');
        xlim([-0.1 0.4]);
        ylim(first_plot_ylim);
        plot_flicker(str2double(freq_interest),temp_ylim,10,condition_color(examples{ex,3}),1);
        ylim(first_plot_ylim2);
        box off;
        xline(0,'--','LineWidth',1,'Label','','LabelHorizontalAlignment','left');
        set(gca,'visible','off');

        ax=gca;
        %ax.Position=[0.1 0.01 0.89 0.98];
        ax.Children(2).FaceColor=[0.5 0.5 0.5];

        %a2=annotation('line',[ax.Position(1) ax.Position(1)+0.1*(ax.Position(3)/range(ax.XLim))],[ax.Position(2) ax.Position(2)],'Linewidth',1);
    
%         temp=num2str(range(ax.YLim));
%         temp2=regexp(temp,'\.');
%         a2=floor(range(ax.YLim)/(10^(temp2-2)))*(10^(temp2-2));
        %a2=annotation('line',[ax.Position(1) ax.Position(1)],[ax.Position(2) ax.Position(2)+temp*(ax.Position(4)/range(ax.YLim))],'Linewidth',1);
        
        
        figure_making('operation','save','filename',[figures_dir '/Figure_S4A_1_LFP.pdf']);
    end

    figure_making('width',6.5,'height',1.5);
    tiledlayout(1,2,'Padding','compact');
    nexttile;
    plot_PSD(examples{ex,3},examples{ex,2},sessions(strcmp({sessions.subjectID},examples{ex,1})).psd_data,sessions(strcmp({sessions.subjectID},examples{ex,1})).psd_data.label,sessions(strcmp({sessions.subjectID},examples{ex,1})).psd_data.condition,condition_color(examples{ex,3}),1,0,1); %plot periodic condition
    hold on;
    plot_PSD('Baseline',examples{ex,2},sessions(strcmp({sessions.subjectID},examples{ex,1})).psd_data,sessions(strcmp({sessions.subjectID},examples{ex,1})).psd_data.label,sessions(strcmp({sessions.subjectID},examples{ex,1})).psd_data.condition,[0 0 0],1,0,1); %plot periodic condition
    l=xline(str2double(freq_interest),'--','LineWidth',1,'Color',condition_color(examples{ex,3}));
    l.Alpha=0.5;
    l=xline(30,'--','LineWidth',1);
    l.Alpha=0.5;
    l=xline(50,'--','LineWidth',1);
    l.Alpha=0.5;
    p=gca;
    temp_ylim=p.YLim;
    %legend(p.Children([7,4]),{'R-AV',condition});
    ax=gca;
    ax.XTick=[0 30 40 50 100];
    %ax.XAxis.FontSize=9;
    xlabel('Frequency (Hz)');
    %ax.YAxis.FontSize=9;
    ylabel('Power (log_1_0(dB))');
    text(85,ax.YLim(end)-range(ax.YLim)/12,['\color[rgb]{' num2str(condition_color(examples{ex,3})) '}' examples{ex,3}],'FontSize',7);
    box off;

    nexttile;
    plot_PSD(['R-' mod_interest],examples{ex,2},sessions(strcmp({sessions.subjectID},examples{ex,1})).psd_data,sessions(strcmp({sessions.subjectID},examples{ex,1})).psd_data.label,sessions(strcmp({sessions.subjectID},examples{ex,1})).psd_data.condition,condition_color(examples{ex,3}),1,0,1); %plot periodic condition
    hold on;
    plot_PSD('Baseline',examples{ex,2},sessions(strcmp({sessions.subjectID},examples{ex,1})).psd_data,sessions(strcmp({sessions.subjectID},examples{ex,1})).psd_data.label,sessions(strcmp({sessions.subjectID},examples{ex,1})).psd_data.condition,[0 0 0],1,0,1); %plot periodic condition
    l=xline(str2double(freq_interest),'--','LineWidth',1,'Color',condition_color(examples{ex,3}));
    l.Alpha=0.5;
    l=xline(30,'--','LineWidth',1);
    l.Alpha=0.5;
    l=xline(50,'--','LineWidth',1);
    l.Alpha=0.5;
    p=gca;
    temp_ylim=p.YLim;
    %legend(p.Children([7,4]),{'R-AV',condition});
    ax=gca;
    ax.XTick=[0 30 40 50 100];
    %ax.XAxis.FontSize=9;
    xlabel('Frequency (Hz)');
    %ax.YAxis.FontSize=9;
    ylabel('Power (log_1_0(dB))');
    text(85,ax.YLim(end)-range(ax.YLim)/12,['\color[rgb]{' num2str(condition_color(['R-' mod_interest])) '}R-' mod_interest],'FontSize',7);
    box off;
    
    figure_making('operation','save','filename',[figures_dir '/' examples{ex,4} '_PSD.pdf']);

end

%Parts of figure: 4.
%2024/02/26

%% define dirs:
root_dir=define_flicker_root_dir;
source_data=[root_dir '/stg-analyses/NatureComm2024-figures/Source Data']; %where source data for figures is stored
figures_dir=[root_dir '/stg-analyses/NatureComm2024-figures'];

%% Panel A- illustration of flicker response in time and frequency domains:

if plot_Figure_4A
    
    examples=readtable([figures_dir '/NatureComm2024_examples.xlsx']);
    examples=table2array(examples(strcmp(examples.example,'Figure_4A'),:));
    examples=[examples(:,2:end) examples(:,1)];
    
    subjectID=examples{1};
    ses=examples{3};
    ch=examples{4};

    %plot example SPEP pulse and subsequent modeled SPEP pulses:
    spep_data=importdata([root_dir '/stg-preproc/sub-' subjectID '/task-spep/ses-' ses '/LFP/spep_erp/sub-' subjectID '_stg-analysis_task-spep_ses-' ses '_nat-spepERP-refLaplacian.mat'],'data_ref_preproc');
    flicker_data=importdata([root_dir '/stg-preproc/sub-' subjectID '/task-flickerneuro/ses-' ses '/LFP/static_ent/sub-' subjectID '_stg-preproc_task-flickerNeuro_ses-' ses '_nat-trialsdata-refLaplacian-preproc.mat'],'data_ref_preproc');

    spep_trials=importdata([root_dir '/stg-preproc/sub-' subjectID '/task-spep/ses-01/sub-' subjectID '_stg-preproc_task-spep_ses-01_nat-beh.mat'],'trials');
    spep_data=importdata([root_dir '/stg-preproc/sub-' subjectID '/task-spep/ses-' ses '/LFP/spep_erp/sub-' subjectID '_stg-analysis_task-spep_ses-' ses '_nat-spepERP-refLaplacian.mat'],'ERP_results_ref_preproc');
    flicker_trials=importdata([root_dir '/stg-preproc/sub-' subjectID '/task-flickerneuro/ses-' ses '/sub-' subjectID '_stg-preproc_task-flickerNeuro_ses-' ses '_nat-beh.mat'],'trials');
    %flicker_data=importdata([root_dir '/stg-preproc/sub-' subjectID '/task-flickerneuro/ses-' ses '/LFP/static_ent/sub-' subjectID '_stg-preproc_task-flickerNeuro_ses-' ses '_nat-trialsdata-refLaplacian-preproc.mat'],'data_ref_preproc');

    %illustrate in time domain:
    %select the trials that you need:
    stim_condition='40Hz-V';
    control_condition='Baseline';
    %keep only stim and control condition trials:
    temp_flicker_trials=flicker_trials;
    temp_index=ismember(temp_flicker_trials.trials_identities(:,1),find(ismember(temp_flicker_trials.condition_code,{stim_condition,control_condition})));
    temp_flicker_trials.trials_identities=temp_flicker_trials.trials_identities(temp_index,:);
    temp_flicker_trials.clinrecording.trials_timestamps=temp_flicker_trials.clinrecording.trials_timestamps(temp_index,:);
    %keep only the control condition trials that we kept in preprocessed data:
    temp_index=[];
    for i=1:size(temp_flicker_trials.trials_identities,1)
        temp=temp_flicker_trials.clinrecording.trials_timestamps(i,1)>=flicker_data.sampleinfo(:,1) & temp_flicker_trials.clinrecording.trials_timestamps(i,2)<=flicker_data.sampleinfo(:,2);
        if any(temp)
            if (sum(temp)~=1)
                error('error');
            else
                temp_index=[temp_index,i];
            end
        end 
    end
    temp_flicker_trials.trials_identities=temp_flicker_trials.trials_identities(temp_index,:);
    temp_flicker_trials.clinrecording.trials_timestamps=temp_flicker_trials.clinrecording.trials_timestamps(temp_index,:);

    [conditionIndex, data]=converttimestamps3(temp_flicker_trials,'clinrecording',2,1); %convert timestamps to timestamps for 2 cycles %NEED TO WORK ON THIS
    conditionIndex(contains(conditionIndex,'occluded'))=[];

    %prepare trials info matrix for processing in fieldtrip (and select a subset of baseline cycles):
    analysis_trials=[];
    for i=conditionIndex
        temp=data{strcmp(conditionIndex,i)};
        temp(:,3)=0;
        temp(:,4)=find(strcmp(conditionIndex,i));
        analysis_trials(end+1:end+size(temp,1),:)=temp;
    end
    analysis_trials=uint32(analysis_trials);

    for i=1:length(flicker_data.trial)
        cfg=[];
        cfg.trials=i;
        cfg.trl=analysis_trials(analysis_trials(:,1)>=flicker_data.sampleinfo(i,1) & analysis_trials(:,2)<=flicker_data.sampleinfo(i,2),:);
        data_withTrials{i}=ft_redefinetrial(cfg,flicker_data);
    end
    cfg=[];
    data_withTrials_sum=ft_appenddata(cfg,data_withTrials{:});


    figure_making('width',2,'height',1.5);

    ch=examples{1,4};

    plot_ERP(spep_data,'V',ch,0,0,0); %plot SPEP trace
    %xlim([0 0.4]); %zoom in x-wise
    curve = findobj(gca,'Type','line');
    curve(1).LineWidth=1;
    curve(1).YData=smooth(curve(1).YData,spep_trials.clinrecording.sampleRate/100);


    %plot modeled flicker response based on linear superposition hypothesis:
    spep_signal=curve(1).YData;
    spep_signal(1:find(curve(1).XData==0))=0;
    spep_signal=resample(spep_signal,40,1); %upsample signal so can divide by 40
    modeled_signal=zeros(1,(0.25+5+1)*spep_data{1}.cfg.previous.previous.previous.previous.hdr.Fs*40); %initialize 0 signal
    modeled_signal(1:find(curve(1).XData==0)*40)=resample(curve(1).YData(1:find(curve(1).XData==0)),40,1);
    current_index=1;
    for i=1:3/(1/40)
        modeled_signal(current_index:current_index+length(spep_signal)-1)=modeled_signal(current_index:current_index+length(spep_signal)-1)+spep_signal; %add spep signal
        current_index=current_index+0.025*spep_data{1}.cfg.previous.previous.previous.previous.hdr.Fs*40; %increment by 25ms
    end
    modeled_signal=resample(modeled_signal,1,40);
    plot(curve(1).XData,modeled_signal(1:length(curve(1).XData)),'Color','k','LineWidth',1);
    axis tight;
    temp=gca;
    ylim([temp.YLim(1)-(temp.YLim(2)-temp.YLim(1))/10 temp.YLim(2)+(temp.YLim(2)-temp.YLim(1))/10]);

    temp=gca;
    ylim([temp.YLim(1) temp.YLim(2)+0.025*(temp.YLim(2)-temp.YLim(1))]);
    temp=gca;
    yrange=[temp.YLim(2)-((temp.YLim(2)-temp.YLim(1))/20) temp.YLim(2)];
    p=patch([0 0.0125 0.0125 0],[yrange(1) yrange(1) yrange(2) yrange(2)],[0.9100 0.4100 0.1700],'EdgeColor','none'); %set patch for SPEP stim pulse
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    p=patch([0.0125 0.025 0.025 0.0125],[yrange(1) yrange(1) yrange(2) yrange(2)],'k','EdgeColor','none'); %set patch for SPEP stim pulse
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

    %draw imaginary pulses and responses:
    start_end=[0 1/40/2;yrange];
    for i=1:(1.3/(1/40))
        if i~=1
            p=patch([start_end(1,:) start_end(1,end:-1:1)],[start_end(2,1) start_end(2,1) start_end(2,2) start_end(2,2)],[0.9100 0.4100 0.1700],'FaceAlpha',0.2,'EdgeColor','none');
            set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            l=plot(curve(1).XData+start_end(1,1),curve(1).YData,'k','LineWidth',1);
            l.Color(4)=0.2;
            start_end(1,:)=start_end(1,:)+1/40/2;
            p=patch([start_end(1,:) start_end(1,end:-1:1)],[start_end(2,1) start_end(2,1) start_end(2,2) start_end(2,2)],'k','FaceAlpha',0.2,'EdgeColor','none');
            set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            start_end(1,:)=start_end(1,:)+1/40/2;
        else
            start_end(1,:)=start_end(1,:)+1/40/2;
            start_end(1,:)=start_end(1,:)+1/40/2;
        end
    end
    set(gca,'children',flipud(get(gca,'children')));

    ax=gca;
    temp=ax.XTickLabel;
    ax.XTickLabel=cellfun(@(x) num2str(str2num(x)*1000),temp,'UniformOutput',false);
    box off;
    xlim([0 0.6]);
    set(gca,'Visible','off');

    ax=gca;
    ax.Position=[0.05 0.01 0.94 0.98];

    a1=annotation('textarrow',[ax.Position(1)-0.025 ax.Position(1)+0.1],[ax.Position(2)+0.05 ax.Position(2)+0.05],'Linewidth',1);
    a2=annotation('textarrow',[ax.Position(1)-0.025 ax.Position(1)-0.025],[ax.Position(2)+0.05 ax.Position(2)+0.5],'Linewidth',1);

    figure_making('operation','save','filename',[figures_dir '/Figure_4A.pdf']);
end


%% fetch data:

all_subjects=fetch_flicker_subjectIDs(root_dir,'spep');

%fetch MNI data:
mni=fetch_mni_anat(root_dir);

%fetch individual subject data:
subject=fetch_subject_data(root_dir,all_subjects,'anat','spep:spep_amp','flickerneuro:ssep_amp');

summary_results=table('Size',[4 4],'VariableTypes',{'string','double','double','double'},'VariableNames',{'Response','V','AV','A'});
summary_results.Response={'None','Flicker only','SPEP only','Flicker and SPEP'}';
for cond={'V','AV','A'}
    
    %create matrix containing result for each contact:
    coresult_mat=[];
    for sub_nber=1:length(subject) %for each subject
            for j=1:size(subject(sub_nber).flickerneuro_ssep_amp_sig,1) %for each electrode
                temp_contact=subject(sub_nber).flickerneuro_ssep_amp_sig.Properties.RowNames{j};
                if ismember(temp_contact,subject(sub_nber).spep_spep_sig.Properties.RowNames) %check if this referenced contact is also present in the spep results
                    temp_coresult=[any(subject(sub_nber).flickerneuro_ssep_amp_sig{j,endsWith(subject(sub_nber).flickerneuro_ssep_amp_sig.Properties.VariableNames,['-' cond{:}])}<0.05) subject(sub_nber).spep_spep_sig{strcmp(temp_contact,subject(sub_nber).spep_spep_sig.Properties.RowNames),cond}<0.05];
                    if all(temp_coresult==[0 0]) %no response to flicker or single pulse
                        dot_color=[0.5 0.5 0.5]; %grey
                        dot_size=1; %small dot
                    elseif all(temp_coresult==[1 0]) %response to flicker but not to single pulse
                        dot_color=[1 0 0]; %red
                        dot_size=2;
                    elseif all(temp_coresult==[1 1]) %response to flicker and single pulse
                        dot_color=[1 0 1]; %purple
                        dot_size=2;
                    elseif all(temp_coresult==[0 1]) %no response to flicker but response to single pulse
                        dot_color=[0 1 1]; %cyan
                        dot_size=2;
                    end
                    temp_contact=strsplit(temp_contact,'-');
                    coresult_mat(end+1,:)=[subject(sub_nber).anat.fs_coords_mni(strcmp(subject(sub_nber).anat.electrodes_info.labels,temp_contact{1}),:) dot_color dot_size];
                end
            end
    end
    
    summary_results(1,cond)={sum(all((coresult_mat(:,4:6)==[0.5 0.5 0.5])'))/size(coresult_mat,1)*100};
    summary_results(2,cond)={sum(all((coresult_mat(:,4:6)==[1 0 0])'))/size(coresult_mat,1)*100};
    summary_results(3,cond)={sum(all((coresult_mat(:,4:6)==[0 1 1])'))/size(coresult_mat,1)*100};
    summary_results(4,cond)={sum(all((coresult_mat(:,4:6)==[1 0 1])'))/size(coresult_mat,1)*100};

end


%make bar graph summarizing the data:
summary_results
disp(round(table2array(summary_results(2:end,2:end))./sum(table2array(summary_results(2:end,2:end)))*100,1)); %display, out of contacts that showed any sensory response, % that show SPEP only, flicker only, or both (for reporting in paper)


figure_making('width',1.5,'height',1.5);
y=summary_results{:,2:end}';
x=1:3;
b=bar(x,y,'stacked');
axis tight;
b(1).FaceColor=[0.5 0.5 0.5]; %grey
b(2).FaceColor=[1 0 0]; %red
b(3).FaceColor=[0 1 1]; %cyan
b(4).FaceColor=[1 0 1]; %purple
temp=gca;
temp.XTickLabel={['\color[rgb]{' num2str(condition_color('V')) '}V'],['\color[rgb]{' num2str(condition_color('AV')) '}AV'],['\color[rgb]{' num2str(condition_color('A')) '}A']};
xlabel('Modality');
ylabel('Percent of contacts');

figure_making('operation','save','filename',[figures_dir '/Figure_4C.pdf']);
writetable(summary_results,[source_data '/Figure_4C.csv']);


%% scatter plot of amplitude of flicker modulation vs SPEP

examples=readtable([figures_dir '/NatureComm2024_examples.xlsx']);
examples=table2array(examples(startsWith(examples.example,'Figure_4D'),:));
examples=[examples(:,2:end) examples(:,1)];

%create correlation matrix of amplitudes normalized by pt and modality
%(i.e. 0 is lowest amplitude for that patient and modality, 1 is highest
%amplitude for that patient and modality):
[~,repo]=define_flicker_root_dir;
path(path,genpath(repo.fieldtrip)); %code below uses the MATLAB normalize function, as opposed to the FieldTrip normalize function; hence, puting the FieldTrip repo to end of path
normalization_range=[0.001 1];
spep_flicker_label={};
spep_flicker_corr_amp=[]; 
spep_flicker_corr_sig=[];
for subject_nber=1:length(subject) %for each subject
    for cond={'V','AV','A'} %for each modality
        temp_spep_amp=subject(subject_nber).spep_spep_amp_val(table2array(subject(subject_nber).spep_spep_sig(:,cond))<0.05,cond); %get all significant SPEP amplitudes for that condition
        temp_spep_amp{:,:}=normalize(table2array(temp_spep_amp),'range',normalization_range); %normalize (range 0,1) SPEP amplitudes for that subject and modality
        
        flicker_cond_index=endsWith(subject(subject_nber).flickerneuro_ssep_amp_sig.Properties.VariableNames,['Hz-' cond{:}]); %get index of columns corresponding to visual condition
        temp_flicker_sig=double(table2array(subject(subject_nber).flickerneuro_ssep_amp_sig(:,flicker_cond_index))<0.05); %get logical matrix of significant flicker results
        temp_flicker_sig(~temp_flicker_sig(:))=NaN; %convert non-significant results to NaN
        temp_flicker_amp=subject(subject_nber).flickerneuro_ssep_amp_val(:,flicker_cond_index);
        temp_flicker_amp{:,:}=temp_flicker_sig.*table2array(temp_flicker_amp); %convert non-significant amplitudes to NaN
        temp_flicker_amp=temp_flicker_amp(~all(isnan(table2array(temp_flicker_amp))'),:); %only keep channels for which there is at least 1 condition that has significant amplitude
        temp_flicker_cond=cell(1,size(temp_flicker_amp,1)); %create cell array that will keep track of which flicker condition gave highest modulation amplitude, for each channel
        for ch=1:size(temp_flicker_amp,1)
            [~,temp]=max(table2array(temp_flicker_amp(ch,:)));
            temp_flicker_cond{ch}=temp_flicker_amp.Properties.VariableNames{temp};
        end
        temp_flicker_amp=table(max(table2array(temp_flicker_amp)')',temp_flicker_cond','VariableNames',{'amplitude','condition'},'RowNames',temp_flicker_amp.Properties.RowNames'); %keep only maximum amplitudes in table
        temp_flicker_amp{:,1}=normalize(table2array(temp_flicker_amp(:,1)),'range',normalization_range); %normalize (range 0,1) flicker amplitudes for that subject and modality
        
        for ch=1:size(temp_spep_amp,1) %for each channel that has significant SPEP amplitude
            temp_contact=temp_spep_amp.Properties.RowNames{ch}; %get channel label
            if ismember(temp_contact,temp_flicker_amp.Properties.RowNames) %check if this referenced contact is also present in the significant flicker results
                spep_flicker_label{end+1}=[subject(subject_nber).subjectID ';' temp_contact ';' temp_flicker_amp{strcmp(temp_contact,temp_flicker_amp.Properties.RowNames),2}{:} ';' num2str(size(spep_flicker_label,2)+1)];
                spep_flicker_corr_amp(end+1,:)=[temp_spep_amp{ch,:} temp_flicker_amp{strcmp(temp_contact,temp_flicker_amp.Properties.RowNames),1}]; %get maximum amplitude for this modality
                spep_flicker_corr_sig(end+1,:)=[subject(subject_nber).spep_spep_sig{strcmp(temp_contact,subject(subject_nber).spep_spep_sig.Properties.RowNames),cond} min(subject(subject_nber).flickerneuro_ssep_amp_sig{strcmp(temp_contact,subject(subject_nber).flickerneuro_ssep_amp_sig.Properties.RowNames),flicker_cond_index})]; %get min p-value for both SPEP and flicker
            end
        end
    end
end

%plot amplitudes:
figure_making('width',2,'height',2);
s=scatter(log10(spep_flicker_corr_amp(:,1)),log10(spep_flicker_corr_amp(:,2)),'filled','k','MarkerFaceAlpha',0.5,'SizeData',5);
assign_scatter_labels(s,spep_flicker_label);
hold on;
plot(log10(normalization_range),log10(normalization_range),'--k');
xlabel(['Normalized single pulse' newline 'response (log_{10} scale)']);
ylabel(['Normalized flicker' newline 'response (log_{10} scale)']);
temp=gca;
temp.XTick=[-3 0];
temp.XTickLabel={'0','1'};
temp.YTick=[-3 0];
temp.YTickLabel={'0','1'};

%mdl=fitlm(log10(spep_flicker_corr_amp(:,1)),log10(spep_flicker_corr_amp(:,2)));
%disp(['Adjusted R-Squared: ' sprintf('%1.2f',mdl.Rsquared.Adjusted) newline 'Estimated Coefficient: ' sprintf('%1.2f',mdl.Coefficients{2,1})]);

temp_examples=join(string([examples(:,1) examples(:,4) examples(:,5)]),';');
examples_index=[find(startsWith(spep_flicker_label,temp_examples(1))) find(startsWith(spep_flicker_label,temp_examples(2)))];
for i=1:length(examples_index)
    scatter(log10(spep_flicker_corr_amp(examples_index(i),1)),log10(spep_flicker_corr_amp(examples_index(i),2)),'filled','r','MarkerFaceAlpha',1,'SizeData',5);
    scatter(log10(spep_flicker_corr_amp(examples_index(i),1)),log10(spep_flicker_corr_amp(examples_index(i),2)),'r','MarkerEdgeAlpha',0.5,'SizeData',60);
end

%find out how many contacts and subjects were used in this plot:
temp_spep_flicker_label=arrayfun(@(x) strsplit(x{:},';'),spep_flicker_label,'UniformOutput',false);
temp_spep_flicker_label=arrayfun(@(x) strjoin(x{:}(1:2),';'),temp_spep_flicker_label,'UniformOutput',false);
temp_spep_flicker_label=unique(temp_spep_flicker_label);
disp(['Total number of contacts: ' num2str(length(temp_spep_flicker_label))]);
temp_spep_flicker_label=arrayfun(@(x) regexprep(x{:},';.+',''),temp_spep_flicker_label,'UniformOutput',false);
disp(['Total number of subjects: ' num2str(length(unique(temp_spep_flicker_label)))]);

%save figure and associated source data:
figure_making('operation','save','filename',[figures_dir '/Figure_4D_2.pdf']);
writetable(cell2table(num2cell(spep_flicker_corr_amp),'VariableNames',{'norm_SPEP_amp','norm_flicker_amp'}),[source_data '/Figure_4D_2.csv']);

%stats of whether amplitudes follow y=x:
disp('ttest:');
%temp=log10(spep_flicker_corr_amp(:,1))-log10(spep_flicker_corr_amp(:,2));
%[h,p,ci,stats]=ttest(temp)
[h,p,ci,stats]=ttest(log10(spep_flicker_corr_amp(:,1)),log10(spep_flicker_corr_amp(:,2)))



%plot significance values:
figure_making('width',2,'height',2);
scatter(spep_flicker_corr_sig(:,1),spep_flicker_corr_sig(:,2),'filled','k','MarkerFaceAlpha',0.5,'SizeData',5);
hold on;
plot([0 0.05],[0 0.05],'--k');
xlim([0 0.05]);
xlabel(['Single pulse response' newline 'significance value']);
ylabel(['Flicker modulation' newline 'significance value']);

%mdl=fitlm(spep_flicker_corr_sig(:,1),spep_flicker_corr_sig(:,2));
%disp(['Adjusted R-Squared: ' sprintf('%1.2f',mdl.Rsquared.Adjusted) newline 'Estimated Coefficient: ' sprintf('%1.2f',mdl.Coefficients{2,1})]);

for i=1:length(examples_index)
    scatter(spep_flicker_corr_sig(examples_index(i),1),spep_flicker_corr_sig(examples_index(i),2),'filled','r','MarkerFaceAlpha',1,'SizeData',5);
    scatter(spep_flicker_corr_sig(examples_index(i),1),spep_flicker_corr_sig(examples_index(i),2),'r','MarkerEdgeAlpha',0.5,'SizeData',60);
end

%save figure and associated source data:
figure_making('operation','save','filename',[figures_dir '/Figure_4D_3.pdf']);
writetable(cell2table(num2cell(spep_flicker_corr_sig),'VariableNames',{'SPEP_sig','flicker_sig'}),[source_data '/Figure_4D_3.csv']);

%stats of whether significance values follow x=y:
disp('ttest:');
% temp=spep_flicker_corr_sig(:,1)-spep_flicker_corr_sig(:,2);
% [h,p]=ttest(temp)
[h,p,ci,stats]=ttest(spep_flicker_corr_sig(:,1),spep_flicker_corr_sig(:,2))


%% Plot example strong SPEP response, low flicker response (PSD), and evolution of flicker response across trial time (to see if explained by decrease in response):

if plot_Figure_4D_1
    examples=examples(:,[1 4 5]);
    control_condition='Baseline';

    PSD_results{1}=importdata([root_dir '/stg-preproc/sub-' examples{1,1} '/task-flickerneuro/ses-01/LFP/static_ent/sub-' examples{1,1} '_stg-analysis_task-flickerneuro_ses-01_nat-psd-refLaplacian.mat'],'PSD_results_ref_preproc');
    PSD_results{2}=importdata([root_dir '/stg-preproc/sub-' examples{2,1} '/task-flickerneuro/ses-01/LFP/static_ent/sub-' examples{2,1} '_stg-analysis_task-flickerneuro_ses-01_nat-psd-refLaplacian.mat'],'PSD_results_ref_preproc');

    flicker_data{1}=importdata([root_dir '/stg-preproc/sub-' examples{1,1} '/task-flickerneuro/ses-01/LFP/static_ent/sub-' examples{1,1} '_stg-preproc_task-flickerNeuro_ses-01_nat-trialsdata-refLaplacian-preproc.mat'],'data_ref_preproc');
    flicker_data{2}=importdata([root_dir '/stg-preproc/sub-' examples{2,1} '/task-flickerneuro/ses-01/LFP/static_ent/sub-' examples{2,1} '_stg-preproc_task-flickerNeuro_ses-01_nat-trialsdata-refLaplacian-preproc.mat'],'data_ref_preproc');

    spep_data{1}=importdata([root_dir '/stg-preproc/sub-' examples{1,1} '/task-spep/ses-01/LFP/spep_erp/sub-' examples{1,1} '_stg-analysis_task-spep_ses-01_nat-spepERP-refLaplacian.mat'],'data_ref_preproc');
    spep_data{2}=importdata([root_dir '/stg-preproc/sub-' examples{2,1} '/task-spep/ses-01/LFP/spep_erp/sub-' examples{2,1} '_stg-analysis_task-spep_ses-01_nat-spepERP-refLaplacian.mat'],'data_ref_preproc');

    for i=1:size(examples,1)
        PSDresults=PSD_results{i};
        time_series=flicker_data{i};

        %calculate flicker ERP for stim and baseline conditions:
        conditions_to_test={examples{i,3},control_condition};

        %calculate flicker ERP for all trials:
        flicker_ERP=cell(1,length(conditions_to_test));
        for condition=conditions_to_test
            cfg=[];
            cfg.trials=find(time_series.trialinfo==find(strcmp(PSDresults.condition,condition)));
            cfg.keeptrials='yes';
            flicker_ERP{strcmp(conditions_to_test,condition)}=ft_timelockanalysis(cfg,time_series); %uses nearest.m function from FieldTrip (not from MATLAB)
        end

        flicker_ERP{end+1}=['ref_method: Laplacian'];
        flicker_ERP{end+1}=conditions_to_test;

        %plot:
        figure_making('width',3,'height',1);
        tiledlayout(1,2,'TileSpacing','none');
        nexttile(1);
        temp=regexprep(examples{i,3},'.+-','');
        plot_ERP(spep_data{i},temp,examples{i,2},1,0,1);
        hold on;
        plot_ERP(spep_data{i},'occluded_AV',examples{i,2},1,0,1);
        xlim([-0.1 0.5]);
        temp2=gca;
        ylim([min(temp2.Children(4).YData) max(temp2.Children(4).YData)]);
        temp1=get(gca,'YLim');
        plot_flicker(40,temp1,0.0125,condition_color(temp),1);
        ylim([temp1(1) temp1(2)+2]);
        box off;
        ylabel([]);
        set(gca,'Visible','off');
        add_scale_bars(gca,0.1,10,0.01); %add scale bars

        nexttile(2);
        %temp_condition=subject(strcmp(subjectIDs,examples{1,1})).ent_amp(examples{1,2},endsWith(subject(strcmp(subjectIDs,examples{1,1})).ent_amp.Properties.VariableNames,['-' examples{1,3}]));
        line_color=condition_color(temp);
        plot_PSD(examples{i,3},examples{i,2},PSDresults,PSDresults.label,PSDresults.condition,line_color,1,0,1);
        hold on;
        plot_PSD('Baseline',examples{i,2},PSDresults,PSDresults.label,PSDresults.condition,'k',1,0,1);
        temp_freq=str2double(regexprep(examples{i,3},'Hz-.+',''));
        xline(temp_freq,':');
        xlim([temp_freq-10 temp_freq+10]);
        %xlabel('Frequency (Hz)');
        %ylabel('Power (log_{10} scale)');
        ylabel(' ');
        box off;
        %add legend:
        ax=gca;
        text(temp_freq+5,ax.YLim(2)-range(ax.YLim)/6,['\color[rgb]{' num2str(line_color) '}' temp '' newline '\color[rgb]{0 0 0}Control']);

        figure_making('operation','save','filename',[figures_dir '/Figure_4D_1_example-' num2str(i) '.pdf']);
    end
end



%% Comparing real flicker data to simulated flicker data (Figure 4E)

ref_method='Laplacian';
for sub_nber=1:length(subject)
    
    %get z score table and p value table for modeled data:
    subject(sub_nber).flickerneuro_ssep_amp_val_model=readtable([root_dir '/stg-analyses/model-superposition/sub-' subject(sub_nber).subjectID '/ses-01/LFP_zscore_table_ref-' ref_method '.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
    subject(sub_nber).flickerneuro_ssep_amp_sig_model=readtable([root_dir '/stg-analyses/model-superposition/sub-' subject(sub_nber).subjectID '/ses-01/LFP_pvalue_table_ref-' ref_method '.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
  
end
 
 
%aggregate all electrodes from all subjects:
comparison_matrix=cell(1,6);
significance_matrix=cell(1,6);
spep_amp_matrix=cell(1,3);
spep_sig_matrix=cell(1,3);
 
for sub_nber=1:length(all_subjects)
        for ch=subject(sub_nber).flickerneuro_ssep_amp_val.Properties.RowNames'
            if ismember(ch,subject(sub_nber).flickerneuro_ssep_amp_val_model.Properties.RowNames') %if that channel was modeled
                comparison_matrix{1}(end+1,1)=subject(sub_nber).flickerneuro_ssep_amp_val{ch,'40Hz-V'};
                comparison_matrix{2}(end+1,1)=subject(sub_nber).flickerneuro_ssep_amp_val_model{ch,'40Hz-V'};
                significance_matrix{1}(end+1,1)=subject(sub_nber).flickerneuro_ssep_amp_sig{ch,'40Hz-V'};
                significance_matrix{2}(end+1,1)=subject(sub_nber).flickerneuro_ssep_amp_sig_model{ch,'40Hz-V'};
                spep_amp_matrix{1}(end+1,1)=subject(sub_nber).spep_spep_amp_val{ch,'V'};
                spep_sig_matrix{1}(end+1,1)=subject(sub_nber).spep_spep_sig{ch,'V'};
                
                comparison_matrix{3}(end+1,1)=subject(sub_nber).flickerneuro_ssep_amp_val{ch,'40Hz-AV'};
                comparison_matrix{4}(end+1,1)=subject(sub_nber).flickerneuro_ssep_amp_val_model{ch,'40Hz-AV'};
                significance_matrix{3}(end+1,1)=subject(sub_nber).flickerneuro_ssep_amp_sig{ch,'40Hz-AV'};
                significance_matrix{4}(end+1,1)=subject(sub_nber).flickerneuro_ssep_amp_sig_model{ch,'40Hz-AV'};
                spep_amp_matrix{2}(end+1,1)=subject(sub_nber).spep_spep_amp_val{ch,'AV'};
                spep_sig_matrix{2}(end+1,1)=subject(sub_nber).spep_spep_sig{ch,'AV'};
                
                comparison_matrix{5}(end+1,1)=subject(sub_nber).flickerneuro_ssep_amp_val{ch,'40Hz-A'};
                comparison_matrix{6}(end+1,1)=subject(sub_nber).flickerneuro_ssep_amp_val_model{ch,'40Hz-A'};
                significance_matrix{5}(end+1,1)=subject(sub_nber).flickerneuro_ssep_amp_sig{ch,'40Hz-A'};
                significance_matrix{6}(end+1,1)=subject(sub_nber).flickerneuro_ssep_amp_sig_model{ch,'40Hz-A'};
                spep_amp_matrix{3}(end+1,1)=subject(sub_nber).spep_spep_amp_val{ch,'A'};
                spep_sig_matrix{3}(end+1,1)=subject(sub_nber).spep_spep_sig{ch,'A'};
            end
        end
end

%only keep results from channels for which we obtained significant
%entrainment in the real data:
for i=1:2:6 %for each real dataset
    comparison_matrix{i}=comparison_matrix{i}(significance_matrix{i}<0.05,:);
    comparison_matrix{i+1}=comparison_matrix{i+1}(significance_matrix{i}<0.05,:);
    spep_amp_matrix{(i+1)/2}=spep_amp_matrix{(i+1)/2}(significance_matrix{i}<0.05,:);
    
    significance_matrix{i}=significance_matrix{i}(significance_matrix{i}<0.05,:);
    significance_matrix{i+1}=significance_matrix{i+1}(significance_matrix{i}<0.05,:);
    spep_sig_matrix{(i+1)/2}=spep_sig_matrix{(i+1)/2}(significance_matrix{i}<0.05,:);
end

%cap fold-change in power at 10:
threshold=10;
comparison_matrix_thresh=comparison_matrix;
for i=1:length(comparison_matrix)
    comparison_matrix_thresh{i}(comparison_matrix_thresh{i}>threshold)=threshold;
end


%display some numbers for text:
disp('Percent sig contacts (simulated over real) for V, AV and V:')
for i=[2 4 6]
    round(sum(significance_matrix{i}<0.05)/length(significance_matrix{i})*100,1)
end

disp('Median values (real V, AV, A:');
for i=[1 3 5]
    round(median(comparison_matrix_thresh{i}),1)
end

disp('Median values (simulated V, AV, A:');
for i=[2 4 6]
    round(median(comparison_matrix_thresh{i}(significance_matrix{i}<0.05)),1)
end

%% draw violin plots:
%CORRECT? (FOR ABOVE AND BELOW)
figure_making('width',3,'height',3);

v=plot_zscore_violin(comparison_matrix_thresh,{'40Hz-V real','40Hz-V model','40Hz-AV real','40Hz-AV model','40Hz-A real','40Hz-A model'},threshold); %uses function violin.m from Violinplot-Matlab
yline(0,'--');
xline(2.5,'-','LineWidth',2);
xline(4.5,'-','LineWidth',2);
plot([1 2],[threshold+1 threshold+1],'k','LineWidth',1);
plot([3 4],[threshold+1 threshold+1],'k','LineWidth',1);
plot([5 6],[threshold+1 threshold+1],'k','LineWidth',1);
axis tight;
ax=gca;
for i=1:length(v)
    v(i).ScatterPlot.SizeData=2;
    if mod(i,2)~=0
        disp('--------------');
        [h,p,ci,stats]=ttest(comparison_matrix{i},comparison_matrix{i+1})
        if p<0.001
            text(i+0.4,threshold+1.5,'***','FontSize',7);
        elseif p<0.01
            text(i+0.4,threshold+1.5,'**','FontSize',7);
        else
            text(i-0.1,threshold+1.5,['* p=' num2str(p)],'FontSize',7);
        end
        text(i-0.1,ax.YLim(1)-0.2,'real','FontSize',7);
        switch i
            case 1
                text(i+0.5,ax.YLim(1)-0.9,['\color[rgb]{' num2str(condition_color('V')) '}V'],'FontSize',7);
            case 3
                text(i+0.5,ax.YLim(1)-0.9,['\color[rgb]{' num2str(condition_color('AV')) '}AV'],'FontSize',7);
            case 5
                text(i+0.5,ax.YLim(1)-0.9,['\color[rgb]{' num2str(condition_color('A')) '}A'],'FontSize',7);
        end
    else
        text(i-0.25,ax.YLim(1)-0.2,'sim','FontSize',7);
        v(i).ScatterPlot.CData(significance_matrix{i}>=0.5,:)=repmat([0.5 0.5 0.5],[sum(significance_matrix{i}>=0.5),1]); %make non-significant model results grey
    end
end

xticklabels([]);
ax.YTick=0:threshold;
ax.YTickLabel(end)={['\geq' num2str(threshold)]};
ylabel(['Flicker response (fold-change)']);
xlabel([newline newline 'Condition']);
% legend('off');
% colormap(flipud(autumn));
% c=colorbar;
% temp=get(gca,'Position');
% c.Position=[c.Position(1:2) c.Position(3)/2 c.Position(4)];
% set(gca,'Position',temp);
% c.TickLabels={'>0','1','2','3','4','5','6','7','8','9','\geq10'};
 
disp(['Number of significant flicker contacts plotted/analyzed: ' num2str(sum(arrayfun(@(x) length(x{:}),comparison_matrix_thresh(:,[1 3 5]))))]);

%save figure and associated source data:
figure_making('operation','save','filename',[figures_dir '/Figure_4E_1.pdf']);

comparison_thresh_tbl=table('Size',[0 3],'VariableNames',{'cond','real-vs-sim','flicker_mod_amp'},'VariableTypes',{'string','string','double'});
for i=1:length(comparison_matrix_thresh)
    if i==1 || i==2
        cond='visual';
    elseif i==3 || i==4
        cond='audiovisual';
    elseif i==5 || i==6
        cond='auditory';
    end
        
    if i==1 || i==3 || i==5
        real_vs_sim='real';
    elseif i==2 || i==4 || i==6
        real_vs_sim='sim';
    end
    
    comparison_thresh_tbl=[comparison_thresh_tbl;repmat({cond},length(comparison_matrix_thresh{i}),1),repmat({real_vs_sim},length(comparison_matrix_thresh{i}),1),num2cell(comparison_matrix_thresh{i})];
end
writetable(comparison_thresh_tbl,[source_data '/Figure_4E.csv']);

%% draw correlation scatter plot:
%NEED TO CHANGE SO THAT SHOW CORRELATION COEFFICIENT
%ADD DIFFERENT COLORS FOR MODEL SIGNIFICANT RESULTS VS NOT?
figure_making('width',3,'height',3); %showing all points (i.e. not differentiating between significant and non-significant flicker/model)
%colors=repmat([0 0 0],[length([comparison_matrix_thresh{2};comparison_matrix_thresh{4};comparison_matrix_thresh{6}]) 1]);
%colors([model_significance_matrix{1};model_significance_matrix{2};model_significance_matrix{3}]==-1,:)=repmat([0.5 0.5 0.5],[sum([model_significance_matrix{1};model_significance_matrix{2};model_significance_matrix{3}]==-1) 1]);
scatter([comparison_matrix_thresh{1};comparison_matrix_thresh{3};comparison_matrix_thresh{5}],[comparison_matrix_thresh{2};comparison_matrix_thresh{4};comparison_matrix_thresh{6}],5,'k','filled','MarkerFaceAlpha',0.5);
temp=gca;
xlabel(['Measured flicker modulation' newline '(fold-change in power)']);
ylabel(['Simulated flicker modulation' newline '(fold-change in power)']);
temp.XLim(1)=min([temp.XLim temp.YLim]);
temp.YLim(1)=min([temp.XLim temp.YLim]);
temp.XLim(2)=max([temp.XLim temp.YLim]);
temp.YLim(2)=max([temp.YLim temp.YLim]);
line(temp.XLim(1):temp.XLim(2),temp.YLim(1):temp.YLim(2),'Color','k','LineStyle','--');
temp.XTickLabel(end)={['\geq' num2str(threshold)]};
temp.YTickLabel(end)={['\geq' num2str(threshold)]};

%estimate whether values are different:
[h,p,ci,stats]=ttest([comparison_matrix{1};comparison_matrix{3};comparison_matrix{5}],[comparison_matrix{2};comparison_matrix{4};comparison_matrix{6}])

%save figure and associated source data:
figure_making('operation','save','filename',[figures_dir '/Figure_4E_2.pdf']);
%source data already saved above

%
%2024/02/26

%% define dirs:

root_dir=define_flicker_root_dir;
source_data=[root_dir '/stg-analyses/NatureComm2024-figures/Source Data']; %where source data for figures is stored
figures_dir=[root_dir '/stg-analyses/NatureComm2024-figures'];

metadata_tbl=readtable([root_dir '/FlickerStudyMetadata.xlsx'],'Sheet','Subjects','PreserveVariableNames',1);
metadata_tbl=metadata_tbl(:,{'Subject_ID'});

mni=fetch_mni_anat(root_dir);

%% fetch data:
all_subjects=fetch_flicker_subjectIDs(root_dir,'all');

%fetch anatomical data:
subject=fetch_subject_data(root_dir,all_subjects,'anat');

%fetch flickerneuro data:
all_subjects=fetch_flicker_subjectIDs(root_dir,'flickerneuro');
temp_subject=fetch_subject_data(root_dir,all_subjects','flickerneuro:ssep_amp');
for sub_nber=1:length(temp_subject)
    subject(strcmp({subject.subjectID},temp_subject(sub_nber).subjectID)).flickerneuro_ssep_amp_sig=temp_subject(sub_nber).flickerneuro_ssep_amp_sig;
end

%fetch spep data:
all_subjects=fetch_flicker_subjectIDs(root_dir,'spep');
temp_subject=fetch_subject_data(root_dir,all_subjects','spep:spep_amp');
for sub_nber=1:length(temp_subject)
    subject(strcmp({subject.subjectID},temp_subject(sub_nber).subjectID)).spep_spep_sig=temp_subject(sub_nber).spep_spep_sig;
end

%fetch flickerfreq data:
all_subjects=fetch_flicker_subjectIDs(root_dir,'flickerfreq');
temp_subject=fetch_subject_data(root_dir,all_subjects','flickerfreq:ssep_amp');
for sub_nber=1:length(temp_subject)
    subject(strcmp({subject.subjectID},temp_subject(sub_nber).subjectID)).flickerfreq_ssep_amp_sig=temp_subject(sub_nber).flickerfreq_ssep_amp_sig;
end

%only keep 1st session for flickerfreq subjects who have multiple sessions (so have only 1 set of channel
%locations to represent):
for sub_nber=1:length(subject)
    if iscell(subject(sub_nber).flickerfreq_ssep_amp_sig)
        subject(sub_nber).flickerfreq_ssep_amp_sig=subject(sub_nber).flickerfreq_ssep_amp_sig{1};
    end
end

total_ch=[];
for sub_nber=1:length(subject)
    if ~isempty(subject(sub_nber).flickerneuro_ssep_amp_sig) %if subject did flickerneuro sessions, count channels from that session
        total_ch=[total_ch size(subject(sub_nber).flickerneuro_ssep_amp_sig,1)];
    else %if no flickerneuro session, subject must have had a flickerfreq session- count channels from that session
        total_ch=[total_ch size(subject(sub_nber).flickerfreq_ssep_amp_sig,1)];
    end
end
total_ch= sum(total_ch);

disp(['Estimated total number of contacts across tasks: ' num2str(total_ch)]);

%% calculate percentage of channels that are in abnormal tissue or SOZ:
% 
% fnames=[];
% fnames.root_dir=root_dir;
% field_names={'flickerneuro_ssep_amp_sig','spep_spep_sig','flickerfreq_ssep_amp_sig'};
% total_chs=0;
% total_patho_chs=0;
% for sub_nber=1:length(subject)
%     fnames.subjectID=subject(sub_nber).subjectID;
%     channels=fetch_channels(fnames,'patho_channels');
%     found_field=0;
%     field_nber=1;
%     while ~found_field
%         if ~isempty(subject(sub_nber).(field_names{field_nber}))
%             found_field=1;
%             temp1=subject(sub_nber).(field_names{field_nber}).Properties.RowNames;
%             total_chs=total_chs+length(temp1);
%             for i=1:length(temp1)
%                 temp2=split(temp1{i},{'-','/'});
%                 if any(ismember(temp2,channels.label(ismember(channels.feature,{'abnormal','soz'}))))
%                     total_patho_chs=total_patho_chs+1;
%                 end
%             end
%         end
%         field_nber=field_nber+1;
%     end
% end
% 
% disp(['Across all channels from all subjects, ' num2str(round((total_chs-total_patho_chs)/total_chs*100,1)) '% (' num2str(total_patho_chs) '/' num2str(total_chs) ' in abnormal tissue/SOZ) of channels used for analysis were NOT in abnormal tissue or SOZ.']);

%% optional removal of contacts in/near abnormal region or SOZ:
% remove_patho_channels=0;
% if remove_patho_channels
%     fnames=[];
%     fnames.root_dir=root_dir;
%     for sub_nber=1:length(subject)
%         fnames.subjectID=subject(sub_nber).subjectID;
%         channels=fetch_channels(fnames,'patho_channels');
%         for field_name={'flickerneuro_ssep_amp_sig','spep_spep_sig','flickerfreq_ssep_amp_sig'}
%             if ~isempty(subject(sub_nber).(field_name{:}))
%                 temp1=subject(sub_nber).(field_name{:}).Properties.RowNames;
%                 row_nbers_excl=[];
%                 for i=1:length(temp1)
%                     temp2=split(temp1{i},{'-','/'});
%                     if any(ismember(temp2,channels.label(ismember(channels.feature,{'abnormal','soz'}))))
%                         row_nbers_excl=[row_nbers_excl,i];
%                     end
%                 end
%                 subject(sub_nber).(field_name{:})(row_nbers_excl,:)=[];
%             end
%         end
%     end
% end

%% plot coverage for each task:

figure_making('width',3.5,'height',2.5);
subplot_nber=1;
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for i={'flickerneuro_ssep_amp_sig','spep_spep_sig','flickerfreq_ssep_amp_sig'}
    nexttile;
    %organize all coordinates and z-values into 1 n x 4 matrix:
    ent_mat=[];
    ent_mat_labels={};
    num_subjects=0;
    for sub_nber=1:length(subject) %for each subject
        if ~isempty(subject(sub_nber).(i{:})) %if have data for that task
            num_subjects=num_subjects+1;
            temp=regexprep(subject(sub_nber).(i{:}).Properties.RowNames,'-.+','');
            ent_mat=[ent_mat;subject(sub_nber).anat.fs_coords_mni(ismember(subject(sub_nber).anat.electrodes_info.labels,temp),:)]; %negative z scores are in grey
            ent_mat_labels=[ent_mat_labels;strcat([subject(sub_nber).subjectID '; '],subject(sub_nber).anat.electrodes_info.labels(ismember(subject(sub_nber).anat.electrodes_info.labels,temp)))]; %NOT SURE THIS GIVES US THE CORRECT MATCHING OF LABELS TO DATA POINTS
        end
    end
    
    %plot brain in 3D:
    if strcmp(i,'flickerneuro_ssep_amp_sig')
        task_name='Flicker 5.5/40/80Hz';
    elseif strcmp(i,'spep_spep_sig')
        task_name='Single-pulse';
    elseif strcmp(i,'flickerfreq_ssep_amp_sig')
        task_name='Flicker 5.5-80Hz range';
    end
    title(['{\bf' task_name '}' newline num2str(num_subjects) ' subjects' newline num2str(size(ent_mat,1)) ' contacts'],'FontWeight','normal','FontSize',6);
    disp([i{:} ': ' num2str(size(ent_mat,1)) ' contacts.']);
    ft_plot_mesh(mni.anat.pial,'facealpha',0.02,'edgecolor','none');
    ft_plot_mesh(mni.anat.mesh_lh,'facealpha',0.1,'facecolor','cortex','edgecolor','none');
    ft_plot_mesh(mni.anat.mesh_rh,'facealpha',0.1,'facecolor','cortex','edgecolor','none');
    hold on;
    lighting gouraud;
    camlight;
    axis tight;

    s=scatter3(ent_mat(:,1),ent_mat(:,2),ent_mat(:,3),'k','filled','SizeData',1,'MarkerFaceAlpha',0.5); %make no ent channels slightly transparent
    assign_scatter_labels(s,ent_mat_labels);
    
    set(gca,'XTick',[], 'YTick', [],'ZTick',[]);
end

figure_making('operation','save','filename',[figures_dir '/Figure_S1B.svg']);

%% show rotating gifs of coverages:
% 
% for i={'flickerneuro_ssep_amp_sig','spep_spep_sig','flickerfreq_ssep_amp_sig'}
%     figure('Position',[10 10 560*3 420*3]);
%     %organize all coordinates and z-values into 1 n x 4 matrix:
%     ent_mat=[];
%     ent_mat_labels={};
%     num_subjects=0;
%     for sub_nber=1:length(subject) %for each subject
%         if ~isempty(subject(sub_nber).(i{:})) %if have data for that task
%             num_subjects=num_subjects+1;
%             temp=regexprep(subject(sub_nber).(i{:}).Properties.RowNames,'-.+','');
%             ent_mat=[ent_mat;subject(sub_nber).anat.fs_coords_mni(ismember(subject(sub_nber).anat.electrodes_info.labels,temp),:)]; %negative z scores are in grey
%             ent_mat_labels=[ent_mat_labels;strcat([subject(sub_nber).subjectID '; '],subject(sub_nber).anat.electrodes_info.labels(ismember(subject(sub_nber).anat.electrodes_info.labels,temp)))]; %NOT SURE THIS GIVES US THE CORRECT MATCHING OF LABELS TO DATA POINTS
%         end
%     end
%     
%     %plot brain in 3D:
%     if strcmp(i,'flickerneuro_ssep_amp_sig')
%         task_name='Flicker 5.5-40-80Hz';
%     elseif strcmp(i,'spep_spep_sig')
%         task_name='Single-pulse';
%     elseif strcmp(i,'flickerfreq_ssep_amp_sig')
%         task_name='Flicker 5.5-80Hz range';
%     end
%     %title(['{\bf' task_name '}' newline num2str(num_subjects) ' subjects' newline num2str(size(ent_mat,1)) ' contacts'],'FontWeight','normal','FontSize',6);
%     disp([i{:} ': ' num2str(size(ent_mat,1)) ' contacts.']);
%     ft_plot_mesh(mni.anat.pial,'facealpha',0.05,'edgecolor','none');
%     ft_plot_mesh(mni.anat.mesh_lh,'facealpha',0.15,'facecolor','cortex','edgecolor','none');
%     ft_plot_mesh(mni.anat.mesh_rh,'facealpha',0.15,'facecolor','cortex','edgecolor','none');
%     hold on;
%     lighting gouraud;
%     camlight;
%     axis tight;
% 
%     s=scatter3(ent_mat(:,1),ent_mat(:,2),ent_mat(:,3),'k','filled','SizeData',30,'MarkerFaceAlpha',0.5); %make no ent channels slightly transparent
%     assign_scatter_labels(s,ent_mat_labels);
%     
%     set(gca,'XTick',[], 'YTick', [],'ZTick',[]);
%     
%     create_gif(gcf,20,180:-1:-180,[figures_dir '/coverage-by-task-' task_name '_3D.gif']);
%     
%     close(gcf);
% end


%% Show whether we have flicker noise, for flickerneuro task.

% fetch flickerneuro data:
all_subjects=fetch_flicker_subjectIDs(root_dir,'flickerneuro');
subject=fetch_subject_data(root_dir,all_subjects','anat','flickerneuro:ssep_amp,psd_data');

%calculate significance and amplitude of ssep for each occluded condition:
for sub_nber=1:length(subject)
    PSD_results=subject(sub_nber).flickerneuro_psd_data; %fetch PSD results for that subject/session

    conditions_of_interest=PSD_results.condition(contains(PSD_results.condition,'occluded')); %find which conditions were occluded
    if ~isempty(conditions_of_interest) %if we do have occluded condition(s) for this subject/session
        control_condition='Baseline'; %define control condition
        zscore_table=zeros(length(PSD_results.label),length(conditions_of_interest)); %initialize modulation amplitude table
        pvalue_table=zeros(length(PSD_results.label),length(conditions_of_interest)); %initialize modulation significance table
        %IS THE BELOW LINE THE RIGHT WAY TO DO THIS?
        rand_baseline_trials_index=randsample(size(PSD_results.data{1,strcmp(PSD_results.condition,control_condition)}{1},1),size(PSD_results.data{1,strcmp(PSD_results.condition,conditions_of_interest(1))}{1},1)); %pick random sample of x (number of occluded trials) of the 15 picked baseline trials, so number of baseline trials matches number of occluded trials
        for i=1:size(zscore_table,2) %for each condition

            freq_interest=str2num(regexprep(conditions_of_interest{i},'occluded_|Hz.+',''));
            [~,temp]=min(abs(PSD_results.data{1,strcmp(PSD_results.condition,conditions_of_interest{i})}{3}-freq_interest)); %find index of frequency closest to frequency of interest- have to do this because sample rate of EDF file not an integer sometimes (error with Natus)
            freq_interest_index=temp;

            for ch=1:size(zscore_table,1) %for each channel
                stim_values=[];
                baseline_values=[];

                for tr=1:size(PSD_results.data{ch,strcmp(PSD_results.condition,conditions_of_interest{i})}{1},1) %for however many number of trials of given condition there are
                    current_stim_value=PSD_results.data{ch,strcmp(PSD_results.condition,conditions_of_interest{i})}{1}(tr,freq_interest_index);
                    current_baseline_value=PSD_results.data{ch,strcmp(PSD_results.condition,control_condition)}{1}(rand_baseline_trials_index(tr),freq_interest_index);

                    stim_values=[stim_values current_stim_value];
                    baseline_values=[baseline_values current_baseline_value];
                end
                zscore_table(ch,i)=(mean(stim_values)/mean(baseline_values))-1;
                pvalue_table(ch,i)=pval_randomshuffle([stim_values' baseline_values'],10000);
            end
        end
        subject(sub_nber).occluded_amp=array2table(zscore_table,'RowNames',PSD_results.label,'VariableNames',conditions_of_interest); %make matrix into table
        subject(sub_nber).occluded_sig=array2table(pvalue_table,'RowNames',PSD_results.label,'VariableNames',conditions_of_interest); %make matrix into table
    end
end


%create table storing significance and amplitude values for non-occluded
%and occluded conditions, for all channels
comp_matrix_tbl=table({''},{''},{''},0,0,0,0,'VariableNames',{'subject_id','contact','condition','significance','amp','significance_occluded','amp_occluded'});
for sub_nber=1:length(subject) %for each subject
    if ~isempty(subject(sub_nber).occluded_sig) %if we have occluded data for that subject
        occluded_conditions=subject(sub_nber).occluded_sig.Properties.VariableNames; %get name of stims we did occluded
        
        for condition=occluded_conditions %for each occluded condition
            comp_matrix_tbl(end+1:end+size(subject(sub_nber).occluded_sig,1),:)=[repmat({subject(sub_nber).subjectID},size(subject(sub_nber).occluded_sig,1),1) ...
                                                                           subject(sub_nber).occluded_sig.Properties.RowNames ...
                                                                           repmat(condition,size(subject(sub_nber).occluded_sig,1),1) ...
                                                                           table2cell(subject(sub_nber).flickerneuro_ssep_amp_sig(:,regexprep(condition,'occluded_',''))) ...
                                                                           table2cell(subject(sub_nber).flickerneuro_ssep_amp_val(:,regexprep(condition,'occluded_',''))) ...
                                                                           table2cell(subject(sub_nber).occluded_sig(:,condition)) ...
                                                                           table2cell(subject(sub_nber).occluded_amp(:,condition))];
        end
    end
end
comp_matrix_tbl(1,:)=[];


%display some numbers:
temp=comp_matrix_tbl(comp_matrix_tbl.significance_occluded<0.05,:);
disp(['Out of all contacts, identified ' num2str((length(unique(strcat(temp.subject_id,'; ',temp.contact)))/length(unique(strcat(comp_matrix_tbl.subject_id,'; ',comp_matrix_tbl.contact))))*100) '% of contacts with a significant occluded condition']);
temp=comp_matrix_tbl(comp_matrix_tbl.significance<0.05,:);
disp(['Out of scenarios where had significant modulation, ' num2str(sum(temp.significance_occluded<0.05)/size(temp,1)*100) '% also had significant modulation in occluded condition.']);
temp=temp(temp.significance_occluded<0.05,:);
disp(['Out of those scenarios, ' num2str(sum(temp.amp>temp.amp_occluded)/size(temp,1)*100) '% showed dose-response behavior.']);


%plot scatter of occluded mod amplitudes vs non-occluded mod amplitudes
for fig={'all_contacts','just_occluded-sig'}
    if strcmp('all_contacts',fig)
        temp=comp_matrix_tbl(comp_matrix_tbl.significance<0.05,:);
        threshold=7;
    elseif strcmp('just_occluded-sig',fig)
        temp=comp_matrix_tbl(comp_matrix_tbl.significance<0.05 & comp_matrix_tbl.significance_occluded<0.05,:);
        threshold=2;
    end
    temp.amp(temp.amp>threshold)=threshold;
    temp.amp_occluded(temp.amp_occluded>threshold)=threshold;
    figure_making('width',1.75,'height',1.75);
    hold on;
    alpha_value=0.5;
    for i=1:size(temp,1)
        if temp{i,'significance_occluded'}<0.05
            s=scatter(temp{i,'amp'},temp{i,'amp_occluded'},2,condition_color(regexprep(temp{i,'condition'},'occluded_','')),'filled','MarkerFaceAlpha',alpha_value,'LineWidth',0.2,'MarkerEdgeColor','k','MarkerEdgeAlpha',alpha_value);
            assign_scatter_labels(s,strcat(temp{i,'subject_id'},'; ',temp{i,'contact'}));
        else
            s=scatter(temp{i,'amp'},temp{i,'amp_occluded'},2,condition_color(regexprep(temp{i,'condition'},'occluded_','')),'filled','MarkerFaceAlpha',alpha_value);
            assign_scatter_labels(s,strcat(temp{i,'subject_id'},'; ',temp{i,'contact'}));
        end
    end
    axis tight;

    xlabel(['Stimulation condition' newline '(fold-change in power)']);
    ylabel(['Occluded condition' newline '(fold-change in power)']);

    temp_gca=gca;
    temp_gca.YLim=[min([temp_gca.YLim temp_gca.XLim]) max([temp_gca.YLim temp_gca.XLim])];
    temp_gca.XLim=temp_gca.YLim;
    temp_gca.XAxis.TickValues=temp_gca.YAxis.TickValues;
    plot(temp_gca.XLim,temp_gca.YLim,'--k');
    temp_gca.XTickLabel{end}=['\geq' temp_gca.XTickLabel{end}];
    temp_gca.YTickLabel{end}=['\geq' temp_gca.YTickLabel{end}];

    %add legend:
    text(temp_gca.XLim(1)+range(temp_gca.XLim)/20,temp_gca.YLim(2)-range(temp_gca.YLim)/7,['\color[rgb]{' num2str(condition_color('V')) '}V' newline '\color[rgb]{' num2str(condition_color('AV')) '}AV' newline '\color[rgb]{' num2str(condition_color('A')) '}A']);
    
    %save figure and associated source data:
    if strcmp(fig,'all_contacts')
        figure_making('operation','save','filename',[figures_dir '/Figure_S1A.pdf']);
        for i=1:size(temp,1)
            temp.subject_id{i}=metadata_tbl{strcmp(metadata_tbl.Subject_ID,temp.subject_id{i}),'Subject_ID'}{:};
        end
        writetable(temp,[source_data '/Figure_S1A.csv']);
    end
       
end




% %get list of contacts that were significantly entrained in occluded
% %conditions, and plot then on MNI brain:
% temp=comp_matrix_tbl(comp_matrix_tbl.significance<0.05 & comp_matrix_tbl.significance_occluded<0.05,:);
% figure_making('width',4,'height',4);
% hold on;
% ft_plot_mesh(mni.anat.pial,'facealpha',0.02,'edgecolor','none');
% ft_plot_mesh(mni.anat.mesh_lh,'facealpha',0.1,'facecolor','cortex','edgecolor','none');
% ft_plot_mesh(mni.anat.mesh_rh,'facealpha',0.1,'facecolor','cortex','edgecolor','none');
% hold on;
% lighting gouraud;
% camlight;
% axis tight;
% for i=1:size(temp,1)
%     sub_nber=find(strcmp({subject.subjectID},temp.subject_id{i}));
%     temp2=regexprep(temp.contact{i},'-.+','');
%     temp2=subject(sub_nber).anat.fs_coords_mni(strcmp(subject(sub_nber).anat.electrodes_info.labels,temp2),:);
%     condition=temp.condition{i};
%     c=condition_color(condition);
%     scatter3(temp2(1),temp2(2),temp2(3),30,c,'filled','MarkerFaceAlpha',0.5);
% end
% 
% figure_making('operation','save','filename',[figures_dir '/flickernoise_3D-brain_significant-occluded.svg']);
% 
% %plot PSD plots for all scenarios:
% figure('Units','normalized','Position',[0 0 1 1]);
% tiledlayout(floor(sqrt(size(temp,1)))+1,floor(sqrt(size(temp,1)))+1);
% for i=1:size(temp,1)
%     nexttile;
%     sub_nber=find(strcmp({subject.subjectID},temp.subject_id{i}));
%     condition=temp.condition{i};
%     c=condition_color(condition);
%     plot_PSD(regexprep(temp.condition{i},'occluded_',''),temp.contact{i},subject(sub_nber).flickerneuro_psd_data,subject(sub_nber).flickerneuro_psd_data.label,subject(sub_nber).flickerneuro_psd_data.condition,c,1,0,1);
%     hold on;
%     plot_PSD(temp.condition{i},temp.contact{i},subject(sub_nber).flickerneuro_psd_data,subject(sub_nber).flickerneuro_psd_data.label,subject(sub_nber).flickerneuro_psd_data.condition,[0.5 0.5 0.5],1,0,1);
%     plot_PSD('Baseline',temp.contact{i},subject(sub_nber).flickerneuro_psd_data,subject(sub_nber).flickerneuro_psd_data.label,subject(sub_nber).flickerneuro_psd_data.condition,'k',1,0,1);
%     xline(40,'--','LineWidth',1);
%     xlim([30 50]);
%     title([temp.subject_id{i} '; ' temp.contact{i}]);
% end
% 


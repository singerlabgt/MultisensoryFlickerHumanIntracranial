%Parts of figures: 3, S3B.
%2024/02/26

%% define dirs and fetch data:
root_dir=define_flicker_root_dir;
source_data=[root_dir '/stg-analyses/NatureComm2024-figures/Source Data']; %where source data for figures is stored
figures_dir=[root_dir '/stg-analyses/NatureComm2024-figures'];

all_subjects=fetch_flicker_subjectIDs(root_dir,'flickerneuro');
threshold=10;

mni=fetch_mni_anat(root_dir);
subject=fetch_subject_data(root_dir,all_subjects','anat','flickerneuro:ssep_amp');


%% summarize response by modality and frequency:

flicker_amp=[]; %will store all significant flicker amplitudes in this table
for subject_nber=1:length(subject)
    temp_flicker_amp=subject(subject_nber).flickerneuro_ssep_amp_val;
    temp_flicker_amp.Properties.RowNames=strcat(subject(subject_nber).subjectID,';',temp_flicker_amp.Properties.RowNames);
    temp_flicker_sig=double(table2array(subject(subject_nber).flickerneuro_ssep_amp_sig)<0.05);
    temp_flicker_sig(~temp_flicker_sig)=NaN;
    temp_flicker_amp{:,:}=table2array(temp_flicker_amp).*temp_flicker_sig;
    flicker_amp=[flicker_amp;temp_flicker_amp];
end

%make matrix of response to modalities:
modality_matrix=table('Size',[size(flicker_amp,1) 3],'VariableNames',{'V','AV','A'},'VariableTypes',{'logical','logical','logical'});
for i=1:size(flicker_amp,1)
    if any(~isnan(table2array(flicker_amp(i,endsWith(flicker_amp.Properties.VariableNames,'-V'))))') %if any response to visual
        modality_matrix{i,'V'}=1;
    end
    
    if any(~isnan(table2array(flicker_amp(i,endsWith(flicker_amp.Properties.VariableNames,'-AV'))))') %if any response to visual
        modality_matrix{i,'AV'}=1;
    end
    
    if any(~isnan(table2array(flicker_amp(i,endsWith(flicker_amp.Properties.VariableNames,'-A'))))') %if any response to visual
        modality_matrix{i,'A'}=1;
    end
end

disp(['Total number of contacts: ' num2str(size(modality_matrix,1)) '.']);

%plot venn diagram of proportion of contacts responding to specific
%modalities:
figure_making('width',1.5,'height',1.5);
venn([sum(modality_matrix{:,'V'}),sum(modality_matrix{:,'AV'}),sum(modality_matrix{:,'A'})],[sum(all(modality_matrix{:,{'V','AV'}}')),sum(all(modality_matrix{:,{'AV','A'}}')),sum(all(modality_matrix{:,{'A','V'}}')),sum(all(modality_matrix{:,{'V','AV','A'}}'))],'FaceColor',{[0.9100 0.4100 0.1700],'g','b'});
axis tight;
axis off;
%add modality annotations:
annotation('textbox', [0, 0.22, 1, 0], 'string', ['V (' num2str(sum(modality_matrix{:,{'V'}})) ')'],'FontWeight','bold','Color',[0.9100 0.4100 0.1700],'EdgeColor','none');
annotation('textbox', [0.65, 0.17, 1, 0], 'string', ['AV (' num2str(sum(modality_matrix{:,{'AV'}})) ')'],'FontWeight','bold','Color','g','EdgeColor','none');
annotation('textbox', [0, 1, 1, 0], 'string', ['A (' num2str(sum(modality_matrix{:,{'A'}})) ')'],'FontWeight','bold','Color','b','EdgeColor','none');
annotation('textbox', [0.5, 1, 1, 0], 'string', ['No mod (' num2str(sum(all((modality_matrix{:,{'V','AV','A'}}==[0 0 0])'))) ')'],'FontWeight','bold','EdgeColor','none');
%add percentages:
annotation('textbox', [0.16, 0.4, 0, 0], 'string', num2str(sum(all((modality_matrix{:,{'V','AV','A'}}==[1 0 0])'))),'FontWeight','bold');
annotation('textbox', [0.4, 0.36, 0, 0], 'string', num2str(sum(all((modality_matrix{:,{'V','AV','A'}}==[1 1 0])'))),'FontWeight','bold');
annotation('textbox', [0.71, 0.5, 0, 0], 'string', num2str(sum(all((modality_matrix{:,{'V','AV','A'}}==[0 1 0])'))),'FontWeight','bold');
annotation('textbox', [0.35, 0.65, 0, 0], 'string', num2str(sum(all((modality_matrix{:,{'V','AV','A'}}==[1 1 1])'))),'FontWeight','bold');
annotation('textbox', [0.17, 0.6, 0, 0], 'string', num2str(sum(all((modality_matrix{:,{'V','AV','A'}}==[1 0 1])'))),'FontWeight','bold');
annotation('textbox', [0.22, 0.9, 0, 0], 'string', num2str(sum(all((modality_matrix{:,{'V','AV','A'}}==[0 0 1])'))),'FontWeight','bold');
annotation('textbox', [0.58, 0.94, 0, 0], 'string', num2str(sum(all((modality_matrix{:,{'V','AV','A'}}==[0 1 1])'))),'FontWeight','bold');
line([6.5 8.5],[12 16.5],'Color','k');

%save figure and associated source data:
figure_making('operation','save','filename',[figures_dir '/Figure_3A_1.pdf']);
writetable(modality_matrix,[source_data '/Figure_3A_1.csv']);

%perform chi-square test of difference between proportions of channels responding to one modality vs another:
temp_cond=[repmat({'V'},size(modality_matrix,1),1);repmat({'AV'},size(modality_matrix,1),1);repmat({'A'},size(modality_matrix,1),1)];
temp_response=~([modality_matrix.V;modality_matrix.AV;modality_matrix.A]);
[tbl,chi2stat,pval]=crosstab(temp_response,temp_cond)

%make matrix of response to frequencies:
frequency_matrix=table('Size',[size(flicker_amp,1) 3],'VariableNames',{'5.5Hz','40Hz','80Hz'},'VariableTypes',{'logical','logical','logical'});
for i=1:size(flicker_amp,1)
    if any(~isnan(table2array(flicker_amp(i,startsWith(flicker_amp.Properties.VariableNames,'5.5Hz'))))') %if any response to visual
        frequency_matrix{i,'5.5Hz'}=1;
    end
    
    if any(~isnan(table2array(flicker_amp(i,startsWith(flicker_amp.Properties.VariableNames,'40Hz'))))') %if any response to visual
        frequency_matrix{i,'40Hz'}=1;
    end
    
    if any(~isnan(table2array(flicker_amp(i,startsWith(flicker_amp.Properties.VariableNames,'80Hz'))))') %if any response to visual
        frequency_matrix{i,'80Hz'}=1;
    end
end

%plot venn diagram of proportion of contacts responding to specific
%frequencies:
figure_making('width',1.5,'height',1.5);
venn([sum(frequency_matrix{:,'5.5Hz'}),sum(frequency_matrix{:,'40Hz'}),sum(frequency_matrix{:,'80Hz'})],[sum(all(frequency_matrix{:,{'5.5Hz','40Hz'}}')),sum(all(frequency_matrix{:,{'40Hz','80Hz'}}')),sum(all(frequency_matrix{:,{'80Hz','5.5Hz'}}')),sum(all(frequency_matrix{:,{'5.5Hz','40Hz','80Hz'}}'))],'FaceColor',{[0.8 0.8 0.8],[0.6 0.6 0.6],[0.4 0.4 0.4]});
axis tight;
axis off;
%add modality annotations:
annotation('textbox', [0, 0.2, 1, 0], 'string', ['5.5Hz (' num2str(sum(frequency_matrix{:,{'5.5Hz'}})) ')'],'FontWeight','bold','Color',[0.8 0.8 0.8],'EdgeColor','none');
annotation('textbox', [0.62, 0.15, 1, 0], 'string', ['40Hz (' num2str(sum(frequency_matrix{:,{'40Hz'}})) ')'],'FontWeight','bold','Color',[0.6 0.6 0.6],'EdgeColor','none');
annotation('textbox', [0, 1, 1, 0], 'string', ['80Hz ' newline '(' num2str(sum(frequency_matrix{:,{'80Hz'}})) ')'],'FontWeight','bold','Color',[0.4 0.4 0.4],'EdgeColor','none');
annotation('textbox', [0.5, 1, 1, 0], 'string', ['No mod (' num2str(sum(all((frequency_matrix{:,{'5.5Hz','40Hz','80Hz'}}==[0 0 0])'))) ')'],'FontWeight','bold','EdgeColor','none');
%add percentages:
annotation('textbox', [0.21, 0.35, 0, 0], 'string', num2str(sum(all((frequency_matrix{:,{'5.5Hz','40Hz','80Hz'}}==[1 0 0])'))),'FontWeight','bold');
annotation('textbox', [0.4, 0.34, 0, 0], 'string', num2str(sum(all((frequency_matrix{:,{'5.5Hz','40Hz','80Hz'}}==[1 1 0])'))),'FontWeight','bold');
annotation('textbox', [0.65, 0.5, 0, 0], 'string', num2str(sum(all((frequency_matrix{:,{'5.5Hz','40Hz','80Hz'}}==[0 1 0])'))),'FontWeight','bold');
annotation('textbox', [0.35, 0.6, 0, 0], 'string', num2str(sum(all((frequency_matrix{:,{'5.5Hz','40Hz','80Hz'}}==[1 1 1])'))),'FontWeight','bold');
annotation('textbox', [0.2, 0.6, 0, 0], 'string', num2str(sum(all((frequency_matrix{:,{'5.5Hz','40Hz','80Hz'}}==[1 0 1])'))),'FontWeight','bold');
annotation('textbox', [0.25, 0.85, 0, 0], 'string', num2str(sum(all((frequency_matrix{:,{'5.5Hz','40Hz','80Hz'}}==[0 0 1])'))),'FontWeight','bold');
annotation('textbox', [0.62, 0.93, 0, 0], 'string', num2str(sum(all((frequency_matrix{:,{'5.5Hz','40Hz','80Hz'}}==[0 1 1])'))),'FontWeight','bold');
line([7 11],[10.5 16],'Color','k');

%save figure and associated source data:
figure_making('operation','save','filename',[figures_dir '/Figure_3A_2.pdf']);
writetable(frequency_matrix,[source_data '/Figure_3A_2.csv']);

%perform chi-square test of difference between proportions of channels responding to one frequency vs another:
temp_cond=[repmat({'5.5Hz'},size(frequency_matrix,1),1);repmat({'40Hz'},size(frequency_matrix,1),1);repmat({'80Hz'},size(frequency_matrix,1),1)];
temp_response=~([frequency_matrix{:,'5.5Hz'};frequency_matrix{:,'40Hz'};frequency_matrix{:,'80Hz'}]);
[tbl,chi2stat,pval]=crosstab(temp_response,temp_cond)

%do analysis of optimal stim frequency, for visual-only vs audio-only modalities:
visual_conditions={'5.5Hz-V','40Hz-V','80Hz-V'};
audio_conditions={'5.5Hz-A','40Hz-A','80Hz-A'};
frequencies=[5.5 40 80];
corr_matrix=[];
for subject_nber=1:length(subject) %for each subject
    for ch=1:size(subject(subject_nber).flickerneuro_ssep_amp_sig,1) %for each channel
        if any(table2array(subject(subject_nber).flickerneuro_ssep_amp_sig(ch,visual_conditions))<0.05) && any(table2array(subject(subject_nber).flickerneuro_ssep_amp_sig(ch,audio_conditions))<0.05) %if any of visual and audio conditions show significant modulation
            [~,max_mod_stimfreq_audio]=max(table2array(subject(subject_nber).flickerneuro_ssep_amp_val(ch,audio_conditions)).*(table2array(subject(subject_nber).flickerneuro_ssep_amp_sig(ch,audio_conditions))<0.05));
            [~,max_mod_stimfreq_visual]=max(table2array(subject(subject_nber).flickerneuro_ssep_amp_val(ch,visual_conditions)).*(table2array(subject(subject_nber).flickerneuro_ssep_amp_sig(ch,visual_conditions))<0.05));
            
            corr_matrix(end+1,:)=[frequencies(max_mod_stimfreq_audio) frequencies(max_mod_stimfreq_visual)];
        end
    end
end

corr_matrix=table(corr_matrix(:,1),corr_matrix(:,2),'VariableNames',{'Audio','Visual'});

disp(['Percent of contacts responding to same frequency regardless of modality: ' num2str(round((sum(corr_matrix{:,1}==corr_matrix{:,2})/size(corr_matrix,1))*100,1)) '%']);

%plot heatmap:
h=heatmap(corr_matrix,'Audio','Visual');
h.YDisplayData=h.YDisplayData(end:-1:1);
h=h.ColorDisplayData;
close(gcf);
h=h/sum(sum(h))*100;
figure_making('width',2,'height',1.5);
htmap=imagesc(h);
caxis([min(min(h)) max(max(h))]);
colormap(flip(gray));
colorbar;
ax=gca;
ax.XTickLabel={'5.5Hz','40Hz','80Hz',''};
ax.YTickLabel={'','5.5Hz','','40Hz','','80Hz',''};
xlabel('Top audio stim frequency (Hz)','Color','b');
ylabel('Top visual stim frequency (Hz)','Color',[0.9100 0.4100 0.1700]);
box off;
ax.XAxis.TickLength=[0 0];
ax.YAxis.TickLength=[0 0];

%save figure and associated source data:
figure_making('operation','save','filename',[figures_dir '/Figure_3A_3.pdf']);
writetable(corr_matrix,[source_data '/Figure_3A_3.csv']);

%% show 1 brain for each condition, and create a rotating gif of the brain:
figure_making('width',3.5,'height',4);
elevation=0;
subplot_nber=1;
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
for i={'5.5Hz-V','5.5Hz-AV','5.5Hz-A','40Hz-V','40Hz-AV','40Hz-A','80Hz-V','80Hz-AV','80Hz-A'}
    nexttile(subplot_nber);
    %organize all coordinates and z-values into 1 n x 4 matrix:
    ent_mat=[];
    source_ent_mat=[];
    for sub_nber=1:length(subject) %for each subject
        for j=1:size(subject(sub_nber).flickerneuro_ssep_amp_val,1) %for each electrode
            temp=[];
            if subject(sub_nber).flickerneuro_ssep_amp_sig{j,i}<=0.05 %significant electrode are larger and in color
                temp1=assign_zscore_color(subject(sub_nber).flickerneuro_ssep_amp_val{j,i},threshold);
                temp2=2;
            else %non-significant electrodes are smaller and grey
                temp1=[0.5 0.5 0.5] ;
                temp2=0.25;
            end
            temp_contact=strsplit(subject(sub_nber).flickerneuro_ssep_amp_val.Properties.RowNames{j},'-');
            ent_mat(end+1,:)=[subject(sub_nber).anat.fs_coords_mni(strcmp(subject(sub_nber).anat.electrodes_info.labels,temp_contact{1}),:) temp1 temp2]; %negative z scores are in grey
            source_ent_mat(end+1,:)=[subject(sub_nber).anat.fs_coords_mni(strcmp(subject(sub_nber).anat.electrodes_info.labels,temp_contact{1}),:) subject(sub_nber).flickerneuro_ssep_amp_val{j,i} subject(sub_nber).flickerneuro_ssep_amp_sig{j,i}];
        end
    end
    
    %noent_mat(end+1,:)=[coordinates(strcmp(electrodes_info.labels,z_score_table.Properties.RowNames{j}),:) z_score_table{j,condition}];
    
    %plot brain in 3D:
    %title(char(condition));
    ft_plot_mesh(mni.anat.pial,'facealpha',0.02,'edgecolor','none');
    ft_plot_mesh(mni.anat.mesh_lh,'facealpha',0.1,'facecolor','cortex','edgecolor','none');
    ft_plot_mesh(mni.anat.mesh_rh,'facealpha',0.1,'facecolor','cortex','edgecolor','none');
    hold on;
    %view(180,20);
    lighting gouraud;
    camlight;
    axis tight;
    if ~isempty(ent_mat)
        temp_noent=ent_mat(ent_mat(:,end)==1,:);
        if ~isempty(temp_noent)
            scatter3(temp_noent(:,1),temp_noent(:,2),temp_noent(:,3),temp_noent(:,7),temp_noent(:,4:6),'filled','MarkerFaceAlpha',0.05); %make no ent channels slightly transparent
        end
        temp_ent=ent_mat(ent_mat(:,end)~=1,:);
        if ~isempty(temp_ent)
            scatter3(temp_ent(:,1),temp_ent(:,2),temp_ent(:,3),temp_ent(:,7),temp_ent(:,4:6),'filled','MarkerFaceAlpha',0.6);
        end
    end
    set(gca,'XTick',[], 'YTick', [],'ZTick',[]);
    colormap(flipud(autumn));
    
    %save figure associated source data:
    writetable(array2table(source_ent_mat,'VariableNames',{'mni_x','mni_y','mni_z','mod_amp','mod_sig'}),[source_data '/Figure_3B_' i{:} '.csv']);

    subplot_nber=subplot_nber+1;
end
hp4 = gca;
hp4=hp4.Position;
c=colorbar;
c.TickLabels=arrayfun(@(x) num2str(x),0:2:threshold,'UniformOutput',false);
c.TickLabels(1)={'0'};
c.TickLabels(end)={['\geq' num2str(threshold)]};
h=gcf;
han=axes(h,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
% annotation('textbox', [0.14, 1, 0, 0], 'string', 'Visual','FontWeight','bold');
% annotation('textbox', [0.3, 1, 0, 0], 'string', 'Audio-visual','FontWeight','bold');
% annotation('textbox', [0.65, 1, 0, 0], 'string', 'Audio','FontWeight','bold');
% annotation('textbox', [0, 0.85, 0, 0], 'string', '5.5Hz','FontWeight','bold');
% annotation('textbox', [0, 0.52, 0, 0], 'string', '40Hz','FontWeight','bold');
% annotation('textbox', [0, 0.22, 0, 0], 'string', '80Hz','FontWeight','bold');

%save figure and associated source data:
figure_making('operation','save','filename',[figures_dir '/Figure_3B.svg']);
%create_gif(h,20,180:-1:-180,[figures_dir '/Human-Flicker-Manuscript/Figures/matlab_figs/Spread-modulation/modulation_across-patients_3D.gif']);

disp(['Total number of contacts: ' num2str(size(ent_mat,1))]);

%% show functional network characterization of flicker modulation, using violing plots:

classifying_table=readtable('anat_cat_yeo2011.csv');
[~,temp]=sort(classifying_table.cog_order);
classifying_table=classifying_table(temp,:);

conditions_of_interest={'5.5Hz-AV','40Hz-AV','80Hz-AV'};

for cond=conditions_of_interest
    disp(['Processing ' cond{:} ' -------------------------------------------------']);
    
    comparison_matrix=cell(1,size(classifying_table,1));
    significance_matrix=cell(1,size(classifying_table,1));

    significant_contacts={};
    %significant_contacts=[];
    for subreg=1:size(classifying_table,1) %for each macro_region
        current_labels=classifying_table{subreg,'network'}; %get wmparc labels of interest
        for sub_nber=1:length(subject) %for each subject
            ch_interest=subject(sub_nber).anat.electrodes_info.labels(endsWith(subject(sub_nber).anat.electrodes_info.anatlabels{:,'fs_yeo2011'},current_labels)); %find channels with these labels
            ch_interest_index=find(startsWith(subject(sub_nber).flickerneuro_ssep_amp_val.Properties.RowNames',strcat(ch_interest,'-')));
            temp_sig=[];
            for ch=ch_interest_index
                comparison_matrix{subreg}(end+1,1)=subject(sub_nber).flickerneuro_ssep_amp_val{ch,cond{:}};
                significance_matrix{subreg}(end+1,1)=subject(sub_nber).flickerneuro_ssep_amp_sig{ch,cond{:}};
                temp_sig=[temp_sig subject(sub_nber).flickerneuro_ssep_amp_sig{ch,cond{:}}];
            end
            significant_contacts{subreg,sub_nber}=[num2str(sum(temp_sig<0.05)) '/' num2str(length(temp_sig))];
            %significant_contacts(subreg,sub_nber)=sum(temp_sig<0.05)/length(temp_sig);
        end
    end

    %disp percent of contacts significantly modulated, by subject:
    array2table(significant_contacts','VariableNames',classifying_table.network,'RowNames',all_subjects)

    %save figure associated source data:
    figure_tbl=table('Size',[0 3],'VariableTypes',{'string','double','double'},'VariableNames',{'func_network','mod_sig','mod_amp'});
    for subreg=1:size(classifying_table,1)
        figure_tbl=[figure_tbl;repmat(classifying_table.network(subreg),length(significance_matrix{subreg}),1),num2cell(significance_matrix{subreg}),num2cell(comparison_matrix{subreg})];
        comparison_matrix{subreg}(comparison_matrix{subreg}>10)=10; %max value will be 10 (for figure)
    end

    %only keep significant contacts:
    for i=1:length(comparison_matrix)
        comparison_matrix{i}=comparison_matrix{i}(significance_matrix{i}<0.05);
    end

    if strcmp(cond,'40Hz-AV')
        figure_making('width',7.5,'height',3);
    else
        figure_making('width',7.5,'height',2.5);
    end
    temp=[1 2 4 7 3 5 6]; %sort violins in always the same order of networks (i.e. from highest number of modulated channels to lowest for the 40Hz-AV condition; 5.5Hz-AV and 80Hz-AV conditions follow the same order
    v=plot_zscore_violin(comparison_matrix(temp),classifying_table.network(temp),threshold);
    for i=1:length(v)
        v(i).ScatterPlot.SizeData=5;
        %v(i).ScatterPlot.CData(significance_matrix{temp(i)}>=0.05,:)=repmat([0.5 0.5 0.5],size(v(i).ScatterPlot.CData(significance_matrix{temp(i)}>=0.05,:),1),1);
        %v(i).ScatterPlot.CData(significance_matrix{temp(i)}<0.05,:)=repmat([0 0 0],size(v(i).ScatterPlot.CData(significance_matrix{temp(i)}<0.05,:),1),1);
        text(i-0.5,threshold+0.5,sprintf('%1.1f%% (n=%d)',(sum(significance_matrix{temp(i)}<0.05)/length(significance_matrix{temp(i)}))*100,length(significance_matrix{temp(i)})),'FontSize',7,'FontWeight','bold');
        %ADD INFORMATION ABOUT TOTAL NUMBER OF CONTACTS IN THAT NETWORK (I.E.
        %N=?)
    end
    %v(2).ScatterPlot.CData(model_significance_matrix{1}==-1,:)=repmat([0.5 0.5 0.5],size(v(2).ScatterPlot.CData(model_significance_matrix{1}==-1,:),1),1);
    %v(4).ScatterPlot.CData(model_significance_matrix{2}==-1,:)=repmat([0.5 0.5 0.5],size(v(4).ScatterPlot.CData(model_significance_matrix{2}==-1,:),1),1);
    %v(6).ScatterPlot.CData(model_significance_matrix{3}==-1,:)=repmat([0.5 0.5 0.5],size(v(6).ScatterPlot.CData(model_significance_matrix{3}==-1,:),1),1);
    yline(0,'--');
    ax=gca;
    %xticklabels([]);
    % for i=1:size(classifying_table,1)
    %     text(i-0.25,ax.YLim(1)-0.75,regexprep(classifying_table.network{i},' ','\n'),'FontSize',7);
    % end
    ylabel(['Flicker response (fold-change in power)']);
    xlabel([newline 'Functional network']);
    legend('off');
    temp=gca;
    temp.YTickLabel(end)={['\geq' num2str(threshold)]};
    
    %save figure and associated source data:
    %save figure and associated source data:
    if strcmp(cond,'40Hz-AV')
        output_filename='Figure_3C';
    elseif strcmp(cond,'5.5Hz-AV')
        output_filename='Figure_S3B_1';
    elseif strcmp(cond,'80Hz-AV')
        output_filename='Figure_S3B_2';
    end
    figure_making('operation','save','filename',[figures_dir '/' output_filename '.pdf']);
    writetable(figure_tbl,[source_data '/' output_filename '.csv']);
    
end

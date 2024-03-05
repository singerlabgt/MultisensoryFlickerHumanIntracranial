% Parts of figures: 1D, 1F, S2, S3A.
%2024/02/26

%% define dirs and fetch data:
root_dir=define_flicker_root_dir;
source_data=[root_dir '/stg-analyses/NatureComm2024-figures/Source Data']; %where source data for figures is stored
figures_dir=[root_dir '/stg-analyses/NatureComm2024-figures']; %figures dir

if ~exist(source_data,'dir')
    mkdir(source_data);
end

%fetch all data we need:
mni=fetch_mni_anat(define_flicker_root_dir);

%fetch electrode localization info and z-score values for all subjects:
all_subjects=fetch_flicker_subjectIDs(root_dir,'flickerneuro');
%subject=fetch_subject_data(root_dir,all_subjects','anat','flickerneuro:ssep_amp,psd_data');
subject=fetch_subject_data(root_dir,all_subjects','anat','flickerneuro:ssep_amp');

examples=readtable([figures_dir '/NatureComm2024_examples.xlsx']);
examples=table2array(examples(startsWith(examples.example,{'Figure_1','Figure_S2','Figure_S3A'}),:));
examples=[examples(:,2:end) examples(:,1)];


%% figure on flicker response in sensory regions:

% determine electrodes of interest (i.e. ones in the visual and auditory regions):
classifying_table=readtable('anat_cat.csv');
visual_regions=classifying_table.subregion_wmparc(strcmp(classifying_table.cog,'visual'));
audio_regions=classifying_table.subregion_wmparc(strcmp(classifying_table.cog,'auditory'));

conditions_of_interest={'40Hz-V','40Hz-A','5.5Hz-V','5.5Hz-A','80Hz-V','80Hz-A'};

top_modulated=table('Size',[10 12],'VariableTypes',repmat({'string','double'},1,6),'VariableNames',reshape([conditions_of_interest;strcat(conditions_of_interest,'_amp')],1,[])); 
for cond=conditions_of_interest
    [visual_mat,audio_mat,visual_mat_labels,audio_mat_labels]=fetch_anat_response(subject,visual_regions,audio_regions,cond{:}); %fetch contact mni coordinates, modulation amplitudes and modulation significance
    
    %present list of top 10 modulated
    if endsWith(cond,'-V')
        temp_mat=visual_mat;
        temp_labels=visual_mat_labels;
    elseif endsWith(cond,'-A')
        temp_mat=audio_mat;
        temp_labels=audio_mat_labels;
    end
    temp_labels(temp_mat(:,5)>=0.05)=[];
    temp_mat(temp_mat(:,5)>=0.05,:)=[];
    [~,temp_index]=sort(temp_mat(:,4),'descend');
    top_modulated.(cond{:})=temp_labels(temp_index(1:10))';
    top_modulated.([cond{:} '_amp'])=temp_mat(temp_index(1:10),4);
    
    
    %plot 3D representation of response:
    threshold=10;
    
    %some numbers for manuscript:
    sig_visual_mat=visual_mat(visual_mat(:,5)<0.05,:);
    sig_audio_mat=audio_mat(audio_mat(:,5)<0.05,:);
    disp(cond{:});
    disp(['Quantity sig contacts in visual areas: ' num2str(size(sig_visual_mat,1)) '/' num2str(size(visual_mat,1)) '(' num2str(round(size(sig_visual_mat,1)/size(visual_mat,1)*100,1)) '%)']);
    disp(['Fold-change- median: ' num2str(round(median(sig_visual_mat(:,4)),1)) '; 25th percentile: ' num2str(round(prctile(sig_visual_mat(:,4),25),1)) '; 75th percentile: ' num2str(round(prctile(sig_visual_mat(:,4),75),1))]);
    disp(['Quantity sig contacts in auditory areas: ' num2str(size(sig_audio_mat,1)) '/' num2str(size(audio_mat,1)) '(' num2str(round(size(sig_audio_mat,1)/size(audio_mat,1)*100,1)) '%)']);
    disp(['Fold-change- median: ' num2str(round(median(sig_audio_mat(:,4)),1)) '; 25th percentile: ' num2str(round(prctile(sig_audio_mat(:,4),25),1)) '; 75th percentile: ' num2str(round(prctile(sig_audio_mat(:,4),75),1))]);
    
    mat=[visual_mat;audio_mat];
    mat_labels=[visual_mat_labels,audio_mat_labels];
    temp_mat=mat(mat(:,5)<0.05,:);
    disp([num2str(sum(temp_mat(:,4)>threshold)) ' contacts had modulation greater than threshold, in ' cond{:} ' condition.']);
    
    figure_making('width',2,'height',2);
    colors=assign_zscore_color(mat(:,4),10);
    ft_plot_mesh(mni.anat.pial,'facealpha',0.02,'edgecolor','none');
    ft_plot_mesh(mni.anat.mesh_lh,'facealpha',0.1,'facecolor','cortex','edgecolor','none');
    ft_plot_mesh(mni.anat.mesh_rh,'facealpha',0.1,'facecolor','cortex','edgecolor','none');
    hold on;
    %view(180,20);
    lighting gouraud;
    camlight;
    axis tight;
    %plot significant contacts:
    temp_mod=mat(mat(:,5)<0.05,:);
    temp_colors=colors(mat(:,5)<0.05,:);
    temp_labels=mat_labels(mat(:,5)<0.05);
    s=scatter3(temp_mod(:,1),temp_mod(:,2),temp_mod(:,3),6,temp_colors,'filled','MarkerFaceAlpha',0.8);
    assign_scatter_labels(s,temp_labels);
    %highlight contact of interest:
    if endsWith(cond,'-V')
        temp_example=examples(strcmp(examples(:,6),'visual-region') & strcmp(examples(:,5),cond),:);
    elseif endsWith(cond,'-A')
        temp_example=examples(strcmp(examples(:,6),'auditory-region') & strcmp(examples(:,5),cond),:);
    end
    temp_example=strcmp(temp_labels,strjoin(temp_example([1,4]),';'));
    scatter3(temp_mod(temp_example,1),temp_mod(temp_example,2),temp_mod(temp_example,3),2,[0.2 0.2 0.2],'filled','MarkerFaceAlpha',1);
    scatter3(temp_mod(temp_example,1),temp_mod(temp_example,2),temp_mod(temp_example,3),50,[0.2 0.2 0.2],'MarkerEdgeAlpha',0.8);
%     %highlight with black circle contact of interest:
%     temp=examples(strcmp(examples(:,3),cond),1:2);
%     temp=strjoin({temp{1},regexprep(temp{2},'-.+','')},';');
%     temp=find(strcmp(temp_labels,temp));
%     scatter3(temp_mod(temp,1),temp_mod(temp,2),temp_mod(temp,3),6,'k','MarkerEdgeAlpha',0.8);
    %plot non-significant contacts:
    temp_mod=mat(mat(:,5)>=0.05,:);
    temp_colors=colors(mat(:,5)>=0.05,:);
    temp_labels=mat_labels(mat(:,5)>=0.05);
    s=scatter3(temp_mod(:,1),temp_mod(:,2),temp_mod(:,3),1,[0.5 0.5 0.5],'filled','MarkerFaceAlpha',0.4);
    assign_scatter_labels(s,temp_labels);
    set(gca,'XTick',[], 'YTick', [],'ZTick',[]);
    colormap(flipud(autumn));
    c=colorbar;
    temp=get(gca,'Position');
    c.Position=[c.Position(1:2) c.Position(3)/2 c.Position(4)];
    set(gca,'Position',temp);
    c.TickLabels=arrayfun(@(x) num2str(x),0:threshold/5:threshold,'UniformOutput',false);
    %c.TickLabels(1)={'>0'};
    c.TickLabels(end)={['\geq' num2str(threshold)]};
    c.Label.String=['Flicker modulation' newline '(fold-change in power)'];
    
    %save figure and associated source data:
    if strcmp(cond,'40Hz-V')
        output_filename='Figure_1D_3D-model';
    elseif strcmp(cond,'40Hz-A')
        output_filename='Figure_1F_3D-model';
    elseif strcmp(cond,'5.5Hz-V')
        output_filename='Figure_S2A_1_3D-model';
    elseif strcmp(cond,'80Hz-V')
        output_filename='Figure_S2A_2_3D-model';
    elseif strcmp(cond,'5.5Hz-A')
        output_filename='Figure_S2B_1_3D-model';
    elseif strcmp(cond,'80Hz-A')
        output_filename='Figure_S2B_2_3D-model';
    end
    %determine number of channels and subjects included in plot:
    all_mat_subjects=arrayfun(@(x) strsplit(x{:},';'),mat_labels,'UniformOutput',false);
    all_mat_subjects=arrayfun(@(x) x{:}{1},all_mat_subjects,'UniformOutput',false);
    all_mat_subjects=unique(all_mat_subjects);
    disp([output_filename ' number of channels, subjects: ' num2str(size(mat,1)) ', ' num2str(length(all_mat_subjects))]);
    figure_making('operation','save','filename',[figures_dir '/' output_filename '.svg']);
    writetable(array2table(mat,'VariableNames',{'mni_x','mni_y','mni_z','mod_amp','mod_sig'}),[source_data '/' output_filename '.csv']);
end


%% flicker response in MTL and PFC:

% determine electrodes of interest (i.e. ones in the HPC and PFC):
classifying_table=readtable('anat_cat.csv');

MTL_names=classifying_table.subregion_wmparc(strcmp(classifying_table.MTL_PFC,'MTL'));
PFC_names=classifying_table.subregion_wmparc(strcmp(classifying_table.MTL_PFC,'PFC'));

conditions_of_interest={'40Hz-AV','5.5Hz-AV','80Hz-AV'};

top_modulated=table('Size',[10 12],'VariableTypes',repmat({'string','double'},1,6),'VariableNames',reshape([strcat(conditions_of_interest,'_MTL'),strcat(conditions_of_interest,'_PFC');strcat(conditions_of_interest,'_MTLamp'),strcat(conditions_of_interest,'_PFCamp')],1,[])); 
for cond=conditions_of_interest

    disp(['Plotting for ' cond{:} ' ---------------------------------------------------------']);
    [MTL_mat,PFC_mat,MTL_mat_labels,PFC_mat_labels]=fetch_anat_response(subject,MTL_names,PFC_names,cond{:});
    
    
    if strcmp(cond,'40Hz-AV')
        %display number of contacts localized to HN and FN:
        disp(['Number of contacts in MTL: ' num2str(size(MTL_mat,1))]);
        disp(['Number of contacts in PFC: ' num2str(size(PFC_mat,1))]);
    end

    %select only significant contacts:
    MTL_mat_sig=MTL_mat(MTL_mat(:,5)<0.05,:);
    MTL_mat_sig_labels=MTL_mat_labels(MTL_mat(:,5)<0.05);
    [~,temp_index]=sort(MTL_mat_sig(:,4),'descend');
    if size(MTL_mat_sig,1)<10
        top_modulated.([cond{:} '_MTL'])(1:size(MTL_mat_sig,1))=MTL_mat_sig_labels(temp_index(1:end))';
        top_modulated.([cond{:} '_MTLamp'])(1:size(MTL_mat_sig,1))=MTL_mat_sig(temp_index(1:end),4);
    else
        top_modulated.([cond{:} '_MTL'])=MTL_mat_sig_labels(temp_index(1:10))';
        top_modulated.([cond{:} '_MTLamp'])=MTL_mat_sig(temp_index(1:10),4);
    end
    PFC_mat_sig=PFC_mat(PFC_mat(:,5)<0.05,:);
    PFC_mat_sig_labels=PFC_mat_labels(PFC_mat(:,5)<0.05);
    [~,temp_index]=sort(PFC_mat_sig(:,4),'descend');
    if size(PFC_mat_sig,1)<10
        top_modulated.([cond{:} '_PFC'])(1:size(PFC_mat_sig,1))=PFC_mat_sig_labels(temp_index(1:end))';
        top_modulated.([cond{:} '_PFCamp'])(1:size(PFC_mat_sig,1))=PFC_mat_sig(temp_index(1:end),4);
    else
        top_modulated.([cond{:} '_PFC'])=PFC_mat_sig_labels(temp_index(1:10))';
        top_modulated.([cond{:} '_PFCamp'])=PFC_mat_sig(temp_index(1:10),4);
    end
    disp(['MTL sig contacts fold-change- median: ' num2str(round(median(MTL_mat_sig(:,4)),1)) ', 25th percentile: ' num2str(round(prctile(MTL_mat_sig(:,4),25),1)) ', 75th percentile: ' num2str(round(prctile(MTL_mat_sig(:,4),75),1))]);
    disp(['PFC sig contacts fold-change- median: ' num2str(round(median(PFC_mat_sig(:,4)),1)) ', 25th percentile: ' num2str(round(prctile(PFC_mat_sig(:,4),25),1)) ', 75th percentile: ' num2str(round(prctile(PFC_mat_sig(:,4),75),1))]);

    if strcmp(cond,'40Hz-AV')
        %pick examples to illustrate:
        [~,temp]=max(MTL_mat_sig(:,4));
        MTL_examples=temp;
        MTL_examples_labels=MTL_mat_sig_labels(temp);
        [~,temp]=sort(abs(MTL_mat_sig(:,4)-median(MTL_mat_sig(:,4)))); %need to pick another example close to median because closest to median is not great
        MTL_examples=[MTL_examples;temp(3)];
        MTL_examples_labels=[MTL_examples_labels,MTL_mat_sig_labels(temp(3))];

        [~,temp]=max(PFC_mat_sig(:,4));
        PFC_examples=temp;
        PFC_examples_labels=PFC_mat_sig_labels(temp);
        [~,temp]=sort(abs(PFC_mat_sig(:,4)-median(PFC_mat_sig(:,4)))); %need to pick another example close to median because closest to median is not great
        PFC_examples=[PFC_examples;temp(4)];
        PFC_examples_labels=[PFC_examples_labels,PFC_mat_sig_labels(temp(4))];
        
        temp=[MTL_examples_labels PFC_examples_labels]';
        temp=arrayfun(@(x) strsplit(x{:},';'),temp,'UniformOutput',false);
        for i=1:length(temp)
            temp2=fetch_subject_data(root_dir,temp{i}(1,1),'flickerneuro:psd_data');
            subject(strcmp({subject.subjectID},temp{i}{1,1})).flickerneuro_psd_data=temp2.flickerneuro_psd_data;
            clear temp2;
        end
    end

    %cap at threshold modulation for better visualization:
    threshold=2;
    MTL_mat_sig_thresh=MTL_mat_sig;
    temp=MTL_mat_sig_thresh(:,4)>threshold;
    disp(['Number of HN contacts more than threshold: ' num2str(sum(temp))]);
    MTL_mat_sig_thresh(temp,4)=threshold;
    PFC_mat_sig_thresh=PFC_mat_sig;
    temp=PFC_mat_sig_thresh(:,4)>threshold;
    disp(['Number of FN contacts more than threshold: ' num2str(sum(temp))]);
    PFC_mat_sig_thresh(temp,4)=threshold;


    % panel D left:
    if strcmp(cond,'40Hz-AV')
        figure_making('width',2.8,'height',2.2);
        subplot(2,3,[2 5]);
        p=plot_zscore_violin({MTL_mat_sig_thresh(:,4),PFC_mat_sig_thresh(:,4)},{'MTL','PFC'},threshold);
        scatter(p(1).ScatterPlot.XData(MTL_examples(1)),p(1).ScatterPlot.YData(MTL_examples(1)),'r');
        scatter(p(1).ScatterPlot.XData(MTL_examples(2)),p(1).ScatterPlot.YData(MTL_examples(2)),'r');
        scatter(p(2).ScatterPlot.XData(PFC_examples(1)),p(2).ScatterPlot.YData(PFC_examples(1)),'r');
        scatter(p(2).ScatterPlot.XData(PFC_examples(2)),p(2).ScatterPlot.YData(PFC_examples(2)),'r');
        p(1).MedianPlot.MarkerFaceAlpha=0;
        p(2).MedianPlot.MarkerFaceAlpha=0;
        text(1,threshold+0.2,sprintf('%1.1f%%\n(n=%d)',(size(MTL_mat_sig,1)/size(MTL_mat,1))*100,size(MTL_mat,1)),'FontSize',7,'HorizontalAlignment','center');
        text(2,threshold+0.2,sprintf('%1.1f%%\n(n=%d)',(size(PFC_mat_sig,1)/size(PFC_mat,1))*100,size(PFC_mat,1)),'FontSize',7,'HorizontalAlignment','center');
        yline(0,'--');
        ax=gca;
        ax.XAxis.FontSize=6;
        ax.YTick(end)=threshold;
        ax.YTickLabel(end)={['\geq' num2str(threshold)]};
        ylabel('Flicker modulation (fold-change in power)','FontSize',6);
        xlabel('Region','FontSize',6);
        ax.Position(4)=ax.Position(4)*0.9;
        ax.Position(2)=ax.Position(2)+ax.Position(2)*0.07;


        subplot(2,3,1);
        temp=strsplit(MTL_examples_labels{1},';');
        channel_name=subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.label(startsWith(subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.label,[temp{2} '-']));
        plot_PSD('40Hz-AV',channel_name{:},subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data,subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.label,...
                        subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.condition,condition_color('40Hz-AV'),1,0,1);
        hold on;
        plot_PSD('Baseline',channel_name{:},subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data,subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.label,...
                        subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.condition,'k',1,0,1);
        axis tight;
        xlim([35 45]);
        xline(40,'--','Color',[0.5 0.5 0.5 0.5]);
        box off;
        temp=gca;
        temp.Position(1)=temp.Position(1)*0.5;
        %add legend:
        ax=gca;
        text(40+0.7,ax.YLim(end)-range(ax.YLim)/6,['\color[rgb]{' num2str(condition_color('40Hz-AV')) '}40Hz-AV' newline '\color[rgb]{0 0 0}' 'Baseline'],'FontSize',6);

        subplot(2,3,4);
        temp=strsplit(MTL_examples_labels{2},';');
        channel_name=subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.label(startsWith(subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.label,[temp{2} '-']));
        plot_PSD('40Hz-AV',channel_name{:},subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data,subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.label,...
                        subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.condition,condition_color('40Hz-AV'),1,0,1);
        hold on;
        plot_PSD('Baseline',channel_name{:},subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data,subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.label,...
                        subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.condition,'k',1,0,1);
        axis tight;
        xlim([35 45]);
        xline(40,'--','Color',[0.5 0.5 0.5 0.5]);
        box off;
        temp=gca;
        temp.Position(1)=temp.Position(1)*0.5;

        subplot(2,3,3);
        temp=strsplit(PFC_examples_labels{1},';');
        channel_name=subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.label(startsWith(subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.label,[temp{2} '-']));
        plot_PSD('40Hz-AV',channel_name{:},subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data,subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.label,...
                        subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.condition,condition_color('40Hz-AV'),1,0,1);
        hold on;
        plot_PSD('Baseline',channel_name{:},subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data,subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.label,...
                        subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.condition,'k',1,0,1);
        axis tight;
        xlim([35 45]);
        xline(40,'--','Color',[0.5 0.5 0.5 0.5]);
        box off;
        temp=gca;
        temp.Position(1)=temp.Position(1)+temp.Position(1)*0.1;

        subplot(2,3,6);
        temp=strsplit(PFC_examples_labels{2},';');
        channel_name=subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.label(startsWith(subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.label,[temp{2} '-']));
        plot_PSD('40Hz-AV',channel_name{:},subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data,subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.label,...
                        subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.condition,condition_color('40Hz-AV'),1,0,1);
        hold on;
        plot_PSD('Baseline',channel_name{:},subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data,subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.label,...
                        subject(strcmp({subject.subjectID},temp(1))).flickerneuro_psd_data.condition,'k',1,0,1);
        axis tight;
        xlim([35 45]);
        xline(40,'--','Color',[0.5 0.5 0.5 0.5]);
        box off;
        temp=gca;
        temp.Position(1)=temp.Position(1)+temp.Position(1)*0.1;
        
        %save figure and associated source data:
        figure_making('operation','save','filename',[figures_dir '/Figure_2D.pdf']);
    end


    %plot 3D representation of response:
    figure_making('width',2,'height',1.5);

    mat=[MTL_mat;PFC_mat];
    mat_labels=[MTL_mat_labels,PFC_mat_labels];
    colors=assign_zscore_color(mat(:,4),threshold);
    ft_plot_mesh(mni.anat.pial,'facealpha',0.02,'edgecolor','none');
    ft_plot_mesh(mni.anat.mesh_lh,'facealpha',0.1,'facecolor','cortex','edgecolor','none');
    ft_plot_mesh(mni.anat.mesh_rh,'facealpha',0.1,'facecolor','cortex','edgecolor','none');
    hold on;
    %view(180,20);
    lighting gouraud;
    camlight;
    axis tight;
    %plot significant contacts:
    temp_mod=mat(mat(:,5)<0.05,:);
    temp_colors=colors(mat(:,5)<0.05,:);
    temp_labels=mat_labels(mat(:,5)<0.05);
    s=scatter3(temp_mod(:,1),temp_mod(:,2),temp_mod(:,3),3,temp_colors,'filled','MarkerFaceAlpha',0.8);
    assign_scatter_labels(s,temp_labels);
    %highlight contact of interest:
    temp_examples=examples(ismember(examples(:,6),{'MTL-region','PFC-region'}) & strcmp(examples(:,5),cond),:);
    if ~isempty(temp_examples)
        temp_examples=ismember(temp_labels,strcat(temp_examples(:,1),';',temp_examples(:,4)));
        scatter3(temp_mod(temp_examples,1),temp_mod(temp_examples,2),temp_mod(temp_examples,3),2,[0.2 0.2 0.2],'filled','MarkerFaceAlpha',1);
        scatter3(temp_mod(temp_examples,1),temp_mod(temp_examples,2),temp_mod(temp_examples,3),50,[0.2 0.2 0.2],'MarkerEdgeAlpha',0.8);
    end
    %plot non-significant contacts:
    temp_mod=mat(mat(:,5)>=0.05,:);
    temp_colors=colors(mat(:,5)>=0.05,:);
    temp_labels=mat_labels(mat(:,5)>=0.05);
    s=scatter3(temp_mod(:,1),temp_mod(:,2),temp_mod(:,3),1,[0.5 0.5 0.5],'filled','MarkerFaceAlpha',0.4);
    assign_scatter_labels(s,temp_labels);
    set(gca,'XTick',[], 'YTick', [],'ZTick',[]);
    colormap(flipud(autumn));
    c=colorbar('Location','southoutside');
    temp=get(gca,'Position');
    c.Position=[c.Position(1:2) c.Position(3) c.Position(4)/2];
    set(gca,'Position',temp);
    c.TickLabels=arrayfun(@(x) num2str(x),0:threshold/2:threshold,'UniformOutput',false);
    c.TickLabels(1)={'0'};
    c.TickLabels(end)={['\geq' num2str(threshold)]};
    c.Label.String=['Flicker modulation' newline '(fold-change in power)'];

%     create_gif(gcf,20,180:-1:-180,[figures_dir '/modulation-higher-cog/MTL-PFC_flicker-response-' cond{:} '_significant_3D.svg']);
%     close(gcf);
                
    %save figure and associated source data:
    if strcmp(cond,'40Hz-AV')
        output_filename='Figure_2C_3D-model';
    elseif strcmp(cond,'5.5Hz-AV')
        output_filename='Figure_S3A_1_3D-model';
    elseif strcmp(cond,'80Hz-AV')
        output_filename='Figure_S3A_2_3D-model';
    end
    %determine number of channels and subjects included in plot:
    all_mat_subjects=arrayfun(@(x) strsplit(x{:},';'),mat_labels,'UniformOutput',false);
    all_mat_subjects=arrayfun(@(x) x{:}{1},all_mat_subjects,'UniformOutput',false);
    all_mat_subjects=unique(all_mat_subjects);
    disp([output_filename ' number of channels, subjects: ' num2str(size(mat,1)) ', ' num2str(length(all_mat_subjects))]);
    figure_making('operation','save','filename',[figures_dir '/' output_filename '.svg']);
    writetable([array2table([repmat({'MTL'},[size(MTL_mat,1),1]);repmat({'PFC'},[size(PFC_mat,1),1])],'VariableNames',{'region'}),...
                    array2table([MTL_mat;PFC_mat],'VariableNames',{'mni_x','mni_y','mni_z','mod_amp','mod_sig'})],[source_data '/' output_filename '.csv']);
end

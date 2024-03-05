%Parts of figure: S7.
%2024/02/26

%% define directories:
root_dir=define_flicker_root_dir;
source_data=[root_dir '/stg-analyses/NatureComm2024-figures/Source Data']; %where source data for figures is stored
figures_dir=[root_dir '/stg-analyses/NatureComm2024-figures'];

%% fetch data:
all_subjects=fetch_flicker_subjectIDs(root_dir,'flickerfreq');
mni=fetch_mni_anat(root_dir);

subject=fetch_subject_data(root_dir,all_subjects','anat','flickerfreq:ssep_amp,ssep_plv');

%% plot results of interest (for figure):
%NEED TO PROOFCHECK

%plot on heatmap to relate to anatomical groups:
classifying_table=readtable('anat_cat.csv');
[~,temp]=sort(classifying_table.cog_order);
classifying_table=classifying_table(temp,:);
classifying_table(ismember(classifying_table.subregion_wmparc,{'corpuscallosum','UnsegmentedWhiteMatter','unknown'}),:)=[];

%concatenate tables of significant results:
visual_mod_tbl={};
audio_mod_tbl={};
for sub_nber=1:length(subject)
    if iscell(subject(sub_nber).flickerfreq_ssep_amp_sig) %means there's multiple sessions
        for s=1:length(subject(sub_nber).flickerfreq_ssep_amp_sig)
            temp=table2array(subject(sub_nber).flickerfreq_ssep_amp_sig{s});
            temp(temp>0.05)=NaN;
            temp(temp<0.05)=1;
            temp=table2array(subject(sub_nber).flickerfreq_ssep_amp_val{s}).*temp;
            temp=array2table(temp,'VariableNames',subject(sub_nber).flickerfreq_ssep_amp_sig{s}.Properties.VariableNames,'RowNames',strcat(subject(sub_nber).subjectID,';',subject(sub_nber).flickerfreq_ssep_amp_sig{s}.Properties.RowNames));
            if contains(subject(sub_nber).flickerfreq_ssep_amp_sig{s}.Properties.VariableNames{1},'V')
                visual_mod_tbl=[visual_mod_tbl;temp];
            elseif contains(subject(sub_nber).flickerfreq_ssep_amp_sig{s}.Properties.VariableNames{1},'A')
                audio_mod_tbl=[audio_mod_tbl;temp];
            end
        end
    else
        temp=table2array(subject(sub_nber).flickerfreq_ssep_amp_sig);
        temp(temp>0.05)=NaN;
        temp(temp<0.05)=1;
        temp=table2array(subject(sub_nber).flickerfreq_ssep_amp_val).*temp;
        temp=array2table(temp,'VariableNames',subject(sub_nber).flickerfreq_ssep_amp_sig.Properties.VariableNames,'RowNames',strcat(subject(sub_nber).subjectID,';',subject(sub_nber).flickerfreq_ssep_amp_sig.Properties.RowNames));
        if contains(subject(sub_nber).flickerfreq_ssep_amp_sig.Properties.VariableNames{1},'V')
            visual_mod_tbl=[visual_mod_tbl;temp];
        elseif contains(subject(sub_nber).flickerfreq_ssep_amp_sig.Properties.VariableNames{1},'A')
            audio_mod_tbl=[audio_mod_tbl;temp];
        end
    end
end

%create heatmap of preferred stimulation frequency by brain region:
%NEED TO DOUBLE-CHECK CORRECT
%plot_colormap=flip(autumn(256*2)); %initialize color map:
plot_colormap=cool(256*2);
num_colors=80*2; %determine the number of colors you will need (i.e. Hz 1 to 80, in increments of 0.5Hz)
temp_index=1:round(size(plot_colormap,1)/num_colors):size(plot_colormap,1);
plot_colormap=plot_colormap(temp_index(1:num_colors),:);
plot_colormap=[1 1 1; plot_colormap]; %add white for contacts that showed no frequency preference (i.e. Hz 0).

%here, select based num stim freqs that led to significant modulation,
%pick top stim freq based on modulation amplitude values
for cond={'V','A'}
    if strcmp(cond,'V')
        mod_tbl=visual_mod_tbl;
    elseif strcmp(cond,'A')
        mod_tbl=audio_mod_tbl;
    end
    
    ent_heatmap=classifying_table;
    heatmap_col=0;
    for subjectID=unique(regexprep(mod_tbl.Properties.RowNames,';.+',''))'
        sub_nber=find(strcmp({subject.subjectID},subjectID)); %find subject number
        temp=mod_tbl(startsWith(mod_tbl.Properties.RowNames,subjectID),:); %subset table to contain data only from that subject
        temp.Properties.RowNames=regexprep(temp.Properties.RowNames,'.+;','');
        heatmap_col=heatmap_col+1;
        for subreg=1:size(classifying_table,1) %for each subregion from wmparc
            current_label=classifying_table{subreg,strcmp(classifying_table.Properties.VariableNames,'subregion_wmparc')}; %find name of current wmparc label


            ent_index=find(endsWith(subject(sub_nber).anat.electrodes_info.anatlabels{:,'fs_aparcaseg'},current_label)); %find indices of electrode contacts that are within that labelled subregion
            ent_index=find(startsWith(temp.Properties.RowNames,subject(sub_nber).anat.electrodes_info.labels(ent_index)));
            if length(ent_index)==0 %means no electrode is in this subregion
                ent_heatmap{subreg,size(classifying_table,2)+heatmap_col}={nan(1)}; %this is going to be grey
            elseif length(ent_index)>0 %means at least 1 or more electrodes found in that subregion
                ent_values=temp{ent_index,:};
                if any(sum(~isnan(ent_values)')>6) %if half or more of the stim freqs lead to significant PLV
                    ent_values=ent_values(sum(~isnan(ent_values)')>6,:); %keep only those modulation values for channels where have half or more of stim frequencies significant
                    [~,temp_index]=max(sum(~isnan(ent_values)')); %take channel that has significant modulation to highest number of stim frequencies %WHAT IF THERE IS A TIE BETWEEN 2 CHANNELS OR MORE?
                    ent_values=ent_values(temp_index,:);
                    [~,temp_index]=max(ent_values); %find index of stim frequency with highest modulation
                    top_freq=str2double(regexprep(temp.Properties.VariableNames(temp_index),'Hz.+',''));
                    ent_heatmap{subreg,size(classifying_table,2)+heatmap_col}={top_freq};
                else
                    ent_heatmap{subreg,size(classifying_table,2)+heatmap_col}={-0.1};
                end
            end
        end
    end

    temp=cell2mat(ent_heatmap{:,size(classifying_table,2)+1:end}); %get data to plot
    temp(temp==-0.1)=0; %assign value 0 to contacts that showed no preference

    heatmap_xlabels=cellstr(string(1:heatmap_col));

    figure_making('width',2.5,'height',4);
    ax1=axes;
    h=heatmap(heatmap_xlabels,ent_heatmap.subregion_wmparc,temp,'Colormap',plot_colormap,'CellLabelColor','none','GridVisible','on','MissingDataColor',[0.8 0.8 0.8],'MissingDataLabel','No coverage','FontSize',12);
    h.FontSize=7;
    xlabel('Sessions');
    temp_ax=gca;
    temp_ax.XDisplayLabels=repmat({''},length(temp_ax.XDisplayLabels),1);
    h.ColorLimits=[0 80];
    %h.GridVisible='off';

    %create customized colobar (to replace heatmap colobar):
    if strcmp(cond,'V')
        ax2=axes; %create 2nd axes where correct colorbar will be plotted
        ax2.Visible='off';
        tested_frequencies=regexprep(subject(1).flickerfreq_ssep_amp_sig.Properties.VariableNames,'Hz-.+',''); %get all the frequencies we tested
        tested_frequencies=sort(arrayfun(@(x) str2num(x{:}),tested_frequencies));
        tested_frequencies_index=int16(tested_frequencies*2)+1; %find index of those frequencies in the colorbar
        plot_colormap2=plot_colormap(tested_frequencies_index,:); %get only colors, from colorbar, coresponding to those frequencies
        colormap(ax2,plot_colormap2); %make this new colormap
        c=colorbar(ax2); %show colormap
        temp=0:1/size(plot_colormap2,1):1; %determine where tick marks should fall (for each of all tested frequencies)
        temp(1)=[];
        temp=temp-1/size(plot_colormap2,1)/2;
        temp=temp([1 4 7 10 13 16 19 22 26]); %choose a subset of frequencies to represent as tickmarks
        c.Ticks=temp;
        temp=tested_frequencies([1 4 7 10 13 16 19 22 26])'; %get actual frequencies to put as labels
        temp=strcat(string(temp),'Hz');
        c.TickLabels=temp';
        c.Position(1)=0.7; %reposition colorbar so overlaps heatmap colobar
        c.Position(2)=0.180;
        c.Position(3)=0.13;
        c.Position(4)=0.75;
    else
       %h.ColorbarVisible='off'; 
    end

    axes_prop=gca;
    x_min=0;
    x_max=h.Position(1)+h.Position(3);
    %x_max=axes_prop.Position(1)+axes_prop.Position(3);
    %x_unit=(x_max-x_min)/(length(conditions)*length(all_subjects));

    y_min=axes_prop.Position(2);
    y_max=axes_prop.Position(2)+axes_prop.Position(4);
    y_unit=(y_max-y_min)/size(classifying_table,1);

    categories=unique(classifying_table.cog,'stable');
    categories=categories(end:-1:1);
    num_rows_completed=0;
    for i=1:length(categories)
        num_subregions=sum(strcmp(classifying_table.cog,categories{i}));
        %annotation('textbox',[0 y_min+num_rows_completed*y_unit 5*x_unit num_subregions*y_unit],'String',categories{i},'FontWeight','bold','FontSize',14,'LineWidth',2,'HorizontalAlignment','center','VerticalAlignment','middle');
        annotation('line',[0 x_max],[y_min+num_rows_completed*y_unit-0.1*y_unit y_min+num_rows_completed*y_unit-0.1*y_unit],'LineWidth',0.5);
        num_rows_completed=num_rows_completed+num_subregions;
    end
    annotation('line',[0 x_max],[y_min+num_rows_completed*y_unit-0.1*y_unit y_min+num_rows_completed*y_unit-0.1*y_unit],'LineWidth',0.5);

    if strcmp(cond,'V')
        annotation('rectangle',[0.725 0.05 0.035 0.05],'FaceColor','w');
        annotation('textbox',[0.725+0.03 0.05 0.5 0.05],'String','No modulation','FontSize',6,'EdgeColor','none');
    end
    
    if strcmp(cond,{'V'})
        file_name='Figure_S7A_1';
    elseif strcmp(cond,{'A'})
        file_name='Figure_S7B_1';
    end
    
    disp([file_name ' number of contacts: ' num2str(size(mod_tbl,1))]);
    
    figure_making('operation','save','filename',[figures_dir '/' file_name '.pdf']);
end


%% plot top modulation frequency on 3D brain:
%NEED TO PROOFCHECK

plot_colormap=cool(256*2);
num_colors=80*2; %determine the number of colors you will need (i.e. Hz 1 to 80, in increments of 0.5Hz)
temp_index=1:round(size(plot_colormap,1)/num_colors):size(plot_colormap,1);
plot_colormap=plot_colormap(temp_index(1:num_colors),:);

for cond={'V','A'}
    if strcmp(cond,'V')
        mod_tbl=visual_mod_tbl;
    elseif strcmp(cond,'A')
        mod_tbl=audio_mod_tbl;
    end
    
    %fetch mni coords and preferred stim freq:
    contacts_mat=[];
    for i=1:size(mod_tbl,1)
        sub_nber=find(strcmp({subject.subjectID},regexprep(mod_tbl.Properties.RowNames{i},';.+','')));
        contact_index=regexprep(mod_tbl.Properties.RowNames{i},'.+;','');
        contact_index=regexprep(contact_index,'-.+','');
        contact_index=find(strcmp(subject(sub_nber).anat.electrodes_info.labels,contact_index));
        if sum(~isnan(mod_tbl{i,:}))>=size(mod_tbl,2)/2 %if more than half stim frequencies have signficant entrainment
            [~,temp_index]=max(mod_tbl{i,:});
            top_freq=str2num(regexprep(mod_tbl.Properties.VariableNames{temp_index},'Hz-.+',''));
        else
            top_freq=0;
        end
        contacts_mat(end+1,:)=[subject(sub_nber).anat.electrodes_info.ecoords_mni_world(contact_index,:) top_freq];
    end

    %plot:
    figure_making('width',3.5,'height',4);
    %figure_making('width',8,'height',8);
    ft_plot_mesh(mni.anat.pial,'facealpha',0.02,'edgecolor','none');
    ft_plot_mesh(mni.anat.mesh_lh,'facealpha',0.1,'facecolor','cortex','edgecolor','none');
    ft_plot_mesh(mni.anat.mesh_rh,'facealpha',0.1,'facecolor','cortex','edgecolor','none');
    hold on;
    %view(180,20);
    lighting gouraud;
    camlight;
    axis tight;
    temp_noent=contacts_mat(contacts_mat(:,4)==0,:);
    if ~isempty(temp_noent)
        scatter3(temp_noent(:,1),temp_noent(:,2),temp_noent(:,3),2,[0.5 0.5 0.5],'filled','MarkerFaceAlpha',0.4); %make no ent channels slightly transparent
    end
    temp_ent=contacts_mat(contacts_mat(:,4)~=0,:);
    if ~isempty(temp_ent)
        scatter3(temp_ent(:,1),temp_ent(:,2),temp_ent(:,3),8,temp_ent(:,4),'filled','MarkerFaceAlpha',0.8);
    end
    set(gca,'XTick',[], 'YTick', [],'ZTick',[]);

    colormap(plot_colormap);
    caxis([0.5 80]);
    c=colorbar;
    c.Visible='off';
    
    contacts_tbl=table('Size',[size(contacts_mat,1) 4],'VariableTypes',{'double','double','double','string'},'VariableNames',{'mni_x','mni_y','mni_z','preferred_freq'});
    contacts_tbl(:,:)=[num2cell(contacts_mat(:,1:3)) cellstr(string(contacts_mat(:,4)))];
    for i=1:size(contacts_mat,1)
        if strcmp(contacts_tbl{i,'preferred_freq'},'0')
            contacts_tbl{i,'preferred_freq'}="no sig mod";
        else
            contacts_tbl{i,'preferred_freq'}=strjoin([contacts_tbl{i,'preferred_freq'},"Hz"],'');
        end
    end
        

%     %create customized colobar (to replace heatmap colobar):
%     ax2=axes; %create 2nd axes where correct colorbar will be plotted
%     ax2.Visible='off';
%     tested_frequencies=regexprep(subject(1).flickerfreq_ssep_amp_sig.Properties.VariableNames,'Hz-.+',''); %get all the frequencies we tested
%     tested_frequencies=sort(arrayfun(@(x) str2num(x{:}),tested_frequencies));
%     tested_frequencies_index=int16(tested_frequencies*2); %find index of those frequencies in the colorbar
%     plot_colormap2=plot_colormap(tested_frequencies_index,:); %get only colors, from colorbar, coresponding to those frequencies
%     colormap(ax2,plot_colormap2); %make this new colormap
%     c=colorbar(ax2); %show colormap
%     temp=0:1/size(plot_colormap2,1):1; %determine where tick marks should fall (for each of all tested frequencies)
%     temp(1)=[];
%     temp=temp-1/size(plot_colormap2,1)/2;
%     temp=temp([1 4 7 10 13 16 19 22 26]); %choose a subset of frequencies to represent as tickmarks
%     c.Ticks=temp;
%     temp=tested_frequencies([1 4 7 10 13 16 19 22 26])'; %get actual frequencies to put as labels
%     temp=strcat(string(temp),'Hz');
%     c.TickLabels=temp';
%     c.Position(2)=0.25;
%     c.Position(4)=0.5;
    %create_gif(gcf,20,180:-1:-180,[figures_dir '/top-stim-freq-anat_3D_' cond{:} '.gif']);
    %close(gcf);
    
    %save figure and associated data:
    if strcmp(cond,{'V'})
        file_name='Figure_S7A_2';
    elseif strcmp(cond,{'A'})
        file_name='Figure_S7B_2';
    end
    
    disp([file_name ' number of contacts: ' num2str(size(contacts_mat,1))]);
    
    figure_making('operation','save','filename',[figures_dir '/' file_name '.svg']);
    writetable(contacts_tbl,[source_data '/' file_name '.csv']);
    
end

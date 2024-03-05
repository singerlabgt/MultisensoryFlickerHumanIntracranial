%Part of figure: S6.
%2024/02/26

%% define directories:
root_dir=define_flicker_root_dir;
source_data=[root_dir '/stg-analyses/NatureComm2024-figures/Source Data']; %where source data for figures is stored
figures_dir=[root_dir '/stg-analyses/NatureComm2024-figures'];

%% draw example SPEP response:

%define which examples we'd like to show:
examples=readtable([figures_dir '/NatureComm2024_examples.xlsx']);
examples=table2array(examples(strcmp(examples.example,'Figure_S6A'),[2,5,6]));

%this is an example of auditory SPEP in primary auditory cortex (Heschl's
%gyrus):
subjectID=examples(1);
channel=examples(2);
condition=examples(3);

%load SPEP data:
ERP_results=importdata([root_dir '/stg-preproc/sub-' subjectID{:} '/task-spep/ses-01/LFP/spep_erp/sub-' subjectID{:} '_stg-analysis_task-spep_ses-01_nat-spepERP-refLaplacian.mat'],'ERP_results_ref_preproc');
      
%plot SPEP:
figure_making('width',1.75,'height',1.75);

plot_ERP(ERP_results,condition,channel,1,0,1);
hold on;
plot_ERP(ERP_results,'occluded_AV',channel,1,0,1);

time_to_plot=[-0.05 0.6]; %in s
[~,index1]=min(abs(ERP_results{1,1}.time-time_to_plot(1)));
[~,index2]=min(abs(ERP_results{1,1}.time-time_to_plot(2)));

xlim([ERP_results{1,1}.time(index1) ERP_results{1,1}.time(index2)]);
ylim([-20 55]);
ax=gca;
ax.XTickLabel=arrayfun(@(x) num2str(x*1000),ax.XTick,'UniformOutput',false);
xline(0,'Label','Pulse','LabelHorizontalAlignment','left','FontSize',6);
xlabel('Time (ms)');
ylabel('Electric Potential (uV)');
set(get(gca,'YLabel'),'Rotation',90);
set(gca,'children',flipud(get(gca,'children')));
line_color=condition_color(condition);
temp1=get(gca,'ylim');
x = [0; 0.0125; 0.0125; 0];
y = [temp1(1); temp1(1); temp1(2); temp1(2)];
plot_flicker(40,temp1,0.0125,line_color,1);

ylim([-12 57]);

box off;
set(gca,'Visible','off');

%add scale bars:
add_scale_bars(gca,0.1,20,0.05);

%add legend:
text(0.4,50,['\color[rgb]{' num2str(condition_color(condition)) '}' condition{:} '' newline '\color[rgb]{' num2str(condition_color('occluded_AV')) '}Occluded']);

figure_making('operation','save','filename',[figures_dir '/Figure_S6A.pdf']);


%% 3D plot of SPEP response in function of modality:

all_subjects=fetch_flicker_subjectIDs(root_dir,'spep');

%fetch MNI data:
mni=fetch_mni_anat(root_dir);

%fetch individual subject data:
subject=fetch_subject_data(root_dir,all_subjects,'anat','spep:spep_amp','flickerneuro:ssep_amp');

figure_making('width',4,'height',1.5);
elevation=0;
subplot_nber=1;
tl_ly=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
threshold=20;
for i={'V','AV','A'}
    nexttile;
    num_contacts_greaterThresh=0;
    %organize all coordinates and z-values into 1 n x 4 matrix:
    ent_mat=[];
    ent_tbl=table('Size',[0 5],'VariableTypes',{'double','double','double','double','double'},'VariableNames',{'mni_x','mni_y','mni_z','significance','max_amplitude'});
    for sub_nber=1:length(subject) %for each subject
        for j=1:size(subject(sub_nber).spep_spep_sig,1) %for each electrode
            if subject(sub_nber).spep_spep_sig{j,i}<0.05 %significant electrode are larger and in color
                temp1=assign_zscore_color(subject(sub_nber).spep_spep_amp_val{j,i},threshold);
                temp2=4;
                if subject(sub_nber).spep_spep_amp_val{j,i}>threshold
                    num_contacts_greaterThresh=num_contacts_greaterThresh+1;
                end
            else %non-significant electrodes are smaller and grey
                temp1=[0.5 0.5 0.5];
                temp2=1;
            end
            temp_contact=strsplit(subject(sub_nber).spep_spep_sig.Properties.RowNames{j},'-');
            ent_mat(end+1,:)=[subject(sub_nber).anat.fs_coords_mni(strcmp(subject(sub_nber).anat.electrodes_info.labels,temp_contact{1}),:) temp1 temp2]; %negative z scores are in grey
            ent_tbl{end+1,:}=[subject(sub_nber).anat.fs_coords_mni(strcmp(subject(sub_nber).anat.electrodes_info.labels,temp_contact{1}),:) subject(sub_nber).spep_spep_sig{j,i} subject(sub_nber).spep_spep_amp_val{j,i}];
        end
    end
    
    disp(['There are ' num2str(num_contacts_greaterThresh) ' contacts with higher than threshold response, in the ' i{:} ' condition.']);
    
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
            scatter3(temp_noent(:,1),temp_noent(:,2),temp_noent(:,3),temp_noent(:,end),temp_noent(:,4:6),'filled','MarkerFaceAlpha',0.4); %make no ent channels slightly transparent
        end
        temp_ent=ent_mat(ent_mat(:,end)~=1,:);
        if ~isempty(temp_ent)
            scatter3(temp_ent(:,1),temp_ent(:,2),temp_ent(:,3),temp_ent(:,end),temp_ent(:,4:6),'filled','MarkerFaceAlpha',0.8);
        end
    end
    set(gca,'XTick',[], 'YTick', [],'ZTick',[]);
    colormap(flipud(autumn));
    
    disp(['Total number of contacts included in the spep ' i{:} ' analysis: ' num2str(size(ent_mat,1))]);
    %save figure associated data:
    writetable(ent_tbl,[source_data '/Figure_S6B_' i{:} '.csv']);

    subplot_nber=subplot_nber+1;
end
hp4 = gca;
hp4=hp4.Position;

%create and save colorbar:
create_colorbar=0;
if create_colorbar
    c=colorbar;
    c.Location='southoutside';
    set(gcf,'Units','Normalized');
    set(gcf,'Position',[0 0 1 1]);
    c.TickLabels=arrayfun(@(x) num2str(x),0:threshold/10:threshold,'UniformOutput',false);
    c.TickLabels(1)={'>0'};
    c.TickLabels(end)={['\geq' num2str(threshold)]};
    c.Label.String='Max absolute peak amplitude (\muV)';
    c.FontSize=8;
    saveas(gcf,[figures_dir '/spep-by-modality_across-patients_3D_colorbar.png']);
    close(gcf);
    % no_sig_legend=[c.Position(1) c.Position(2)-0.07 c.Position(3) 0.05];
    % annotation('rectangle',no_sig_legend,'FaceColor',[0.5 0.5 0.5]);
    % no_sig_legend(1)=no_sig_legend(1)+c.Position(3);
    % no_sig_legend(2)=no_sig_legend(2)+c.Position(3);
    % no_sig_legend(3)=1;
    %annotation('textbox',no_sig_legend,'String','No SPEP','EdgeColor','none','FontSize',7);
end

h=gcf;
han=axes(h,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
% annotation('textbox', [0.22, 1, 0, 0], 'string', 'visual','FontWeight','bold');
% annotation('textbox', [0.47, 1, 0, 0], 'string', 'audio-visual','FontWeight','bold');
% annotation('textbox', [0.75, 1, 0, 0], 'string', 'audio','FontWeight','bold');

%save figure:
figure_making('operation','save','filename',[figures_dir '/Figure_S6B.svg']);


%% plot flicker vs spep response, for visual and audio modalities:

%plot 3D representation of flicker vs SPEP response:
figure_making('width',4,'height',1.5);
%figure_making('width',5,'height',2.5);
subplot_nber=1;
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
summary_results=table('Size',[4 4],'VariableTypes',{'string','double','double','double'},'VariableNames',{'Response','V','AV','A'});
summary_results.Response={'None','Flicker only','SPEP only','Flicker and SPEP'}';
for cond={'V','AV','A'}
    
    nexttile;
    
    %create matrix containing result for each contact:
    coresult_mat=[];
    dot_color=[];
    dot_resp='';
    coresult_tbl=table('Size',[0 4],'VariableTypes',{'double','double','double','string',},'VariableNames',{'mni_x','mni_y','mni_z','response'});
    label_mat={};
    for sub_nber=1:length(subject) %for each subject
            for j=1:size(subject(sub_nber).flickerneuro_ssep_amp_sig,1) %for each electrode
                temp_contact=subject(sub_nber).flickerneuro_ssep_amp_sig.Properties.RowNames{j};
                if ismember(temp_contact,subject(sub_nber).spep_spep_sig.Properties.RowNames) %check if this referenced contact is also present in the spep results
                    temp_coresult=[any(subject(sub_nber).flickerneuro_ssep_amp_sig{j,endsWith(subject(sub_nber).flickerneuro_ssep_amp_sig.Properties.VariableNames,['-' cond{:}])}<0.05) subject(sub_nber).spep_spep_sig{strcmp(temp_contact,subject(sub_nber).spep_spep_sig.Properties.RowNames),cond}<0.05];
                    if all(temp_coresult==[0 0]) %no response to flicker or single pulse
                        dot_color=[0.5 0.5 0.5]; %grey
                        dot_size=1; %small dot
                        dot_resp='none';
                    elseif all(temp_coresult==[1 0]) %response to flicker but not to single pulse
                        dot_color=[1 0 0]; %red
                        dot_size=4;
                        dot_resp='flicker';
                    elseif all(temp_coresult==[1 1]) %response to flicker and single pulse
                        dot_color=[1 0 1]; %purple
                        dot_size=4;
                        dot_resp='flicker and single pulse';
                    elseif all(temp_coresult==[0 1]) %no response to flicker but response to single pulse
                        dot_color=[0 1 1]; %cyan
                        dot_size=4;
                        dot_resp='single pulse';
                    end
                    temp_contact=strsplit(temp_contact,'-');
                    coresult_mat(end+1,:)=[subject(sub_nber).anat.fs_coords_mni(strcmp(subject(sub_nber).anat.electrodes_info.labels,temp_contact{1}),:) dot_color dot_size];
                    label_mat{end+1}=[subject(sub_nber).subjectID '; ' temp_contact{1}];
                    coresult_tbl(end+1,:)=[num2cell(subject(sub_nber).anat.fs_coords_mni(strcmp(subject(sub_nber).anat.electrodes_info.labels,temp_contact{1}),:)) dot_resp];
                end
            end
    end
    
    %plot brain in 3D:
    ft_plot_mesh(mni.anat.pial,'facealpha',0.02,'edgecolor','none');
    ft_plot_mesh(mni.anat.mesh_lh,'facealpha',0.1,'facecolor','cortex','edgecolor','none');
    ft_plot_mesh(mni.anat.mesh_rh,'facealpha',0.1,'facecolor','cortex','edgecolor','none');
    hold on;
    %view(180,20);
    lighting gouraud;
    camlight;
    axis tight;
    
    %plot no response contacts:
    temp_mat=coresult_mat(all((coresult_mat(:,4:6)==[0.5 0.5 0.5])'),:);
    temp_labels=label_mat(all((coresult_mat(:,4:6)==[0.5 0.5 0.5])'));
    s=scatter3(temp_mat(:,1),temp_mat(:,2),temp_mat(:,3),temp_mat(:,7),temp_mat(:,4:6),'filled','MarkerFaceAlpha',0.4);
    row = dataTipTextRow('Label',temp_labels);
    s.DataTipTemplate.DataTipRows(end+1) = row;
    
    %plot response contacts:
    temp_mat=coresult_mat(~all((coresult_mat(:,4:6)==[0.5 0.5 0.5])'),:);
    temp_labels=label_mat(~all((coresult_mat(:,4:6)==[0.5 0.5 0.5])'));
    s=scatter3(temp_mat(:,1),temp_mat(:,2),temp_mat(:,3),temp_mat(:,7),temp_mat(:,4:6),'filled','MarkerFaceAlpha',0.6);
    row = dataTipTextRow('Label',temp_labels);
    s.DataTipTemplate.DataTipRows(end+1) = row;
    
    
    summary_results(1,cond)={sum(all((coresult_mat(:,4:6)==[0.5 0.5 0.5])'))/size(coresult_mat,1)*100};
    summary_results(2,cond)={sum(all((coresult_mat(:,4:6)==[1 0 0])'))/size(coresult_mat,1)*100};
    summary_results(3,cond)={sum(all((coresult_mat(:,4:6)==[0 1 1])'))/size(coresult_mat,1)*100};
    summary_results(4,cond)={sum(all((coresult_mat(:,4:6)==[1 0 1])'))/size(coresult_mat,1)*100};
    
    set(gca,'XTick',[], 'YTick', [],'ZTick',[]);
    
    disp(['Total number of contacts included in the ' cond{:} ' analysis: ' num2str(size(coresult_mat,1))]);
    
    %save figure associated data:
    writetable(coresult_tbl,[source_data '/Figure_S6C_' cond{:} '.csv']);
    
end

h=gcf;
han=axes(h,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
% annotation('textbox', [0.2, 1, 0, 0], 'string', 'visual','FontWeight','bold');
% annotation('textbox', [0.45, 1, 0, 0], 'string', 'audio-visual','FontWeight','bold');
% annotation('textbox', [0.75, 1, 0, 0], 'string', 'audio','FontWeight','bold');

% hold on
% temp = gobjects(4,1);
% temp(1) = scatter3(nan,nan,nan,1,[0.5 0.5 0.5],'filled');
% temp(2) = scatter3(nan,nan,nan,8,[1 0 0],'filled');
% temp(3) = scatter3(nan,nan,nan,8,[0 1 1],'filled');
% temp(4) = scatter3(nan,nan,nan,8,[1 0 1],'filled');

% lg=legend(temp,'None','Flicker only','SPEP only','Flicker and SPEP');
% title(lg,'Response');
% lg.Position=[0.83,0.04,0.17,0.19];
% lg.Box='off';

figure_making('operation','save','filename',[figures_dir '/Figure_S6C.svg']);

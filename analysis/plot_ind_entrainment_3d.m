%Plots flicker modulation significance and amplitude on 3D model of subject's brain.
%2024/02/25

function plot_ind_entrainment_3d(fnames,ref_method)
    %initialize subject structure:
    subject=fetch_subject_data(fnames.root_dir,{fnames.subjectID},'anat');
    subject_nber=1;
    
    %get mni template data:
    mni=fetch_mni_anat(fnames.root_dir);
    
    %get flicker modulation significance and amplitude values:
    subject(subject_nber).pvalue=readtable([fnames.root_dir '/stg-analyses/task-' fnames.task '/sub-' fnames.subjectID '/ses-' fnames.ses '/LFP/static_ent/LFP_pvalue_table_ref' ref_method '.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
    subject(subject_nber).zscore=readtable([fnames.root_dir '/stg-analyses/task-' fnames.task '/sub-' fnames.subjectID '/ses-' fnames.ses '/LFP/static_ent/LFP_zscore_table_ref' ref_method '.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
      
    %perpare plotting:
    h=figure('units','normalized','outerposition',[0 0 0.7 1]);
    threshold=10;
    elevation=0;
    %define conditions of interest:
    if strcmp(fnames.task,'flickerneuro') %for flickerneuro task
        conditions={'5.5Hz-V','5.5Hz-AV','5.5Hz-A','40Hz-V','40Hz-AV','40Hz-A','80Hz-V','80Hz-AV','80Hz-A'};
        tiledlayout(3,3,'TileSpacing','compact');
    elseif strcmp(fnames.task,'flickerfreq') %for flickerfreq task
        temp=subject(1).zscore.Properties.VariableNames;
        temp=regexprep(temp,'Hz.+','');
        temp=arrayfun(@(x) str2num(x{:}),temp);
        [~,temp]=sort(temp);
        conditions=subject(1).zscore.Properties.VariableNames(temp);
        tiledlayout(5,6,'TileSpacing','compact');
    end
    
    %plot modulation for each condition:
    subplot_nber=1;
    for i=conditions %for each condition
        %organize electrodes' xyz, modulation amplitude and significance values into an n x 5 matrix:
        ent_mat=[];
        for j=1:size(subject(1).zscore,1) %for each electrode
            temp=[];
            if subject(1).pvalue{j,i}<=0.05 %significant electrode are larger and in color
                temp1=assign_zscore_color(subject(1).zscore{j,i},threshold); %assign color based on modulation amplitude
                temp2=20;
            else %non-significant electrodes are smaller and grey
                temp1=[0.5 0.5 0.5]; %grey
                temp2=1;
            end
            temp_contact=strsplit(subject(1).zscore.Properties.RowNames{j},'-');
            ent_mat(end+1,:)=[subject(1).anat.fs_coords_mni(strcmp(subject(1).anat.electrodes_info.labels,temp_contact{1}),:) temp1 temp2];
        end
        
        %plot modulation on 3D brain:
        nexttile;
        hold on;
        ft_plot_mesh(mni.anat.pial,'facealpha',0.02,'edgecolor','none'); %plot pial surface
        ft_plot_mesh(mni.anat.mesh_lh,'facealpha',0.1,'facecolor','cortex','edgecolor','none'); %plot left hippocampus
        ft_plot_mesh(mni.anat.mesh_rh,'facealpha',0.1,'facecolor','cortex','edgecolor','none'); %plot right hippocampus
        lighting gouraud;
        camlight;
        axis tight;
        if ~isempty(ent_mat)
            %identify non-significant electrodes and plot:
            temp_noent=ent_mat(ent_mat(:,end)==1,:);
            if ~isempty(temp_noent)
                scatter3(temp_noent(:,1),temp_noent(:,2),temp_noent(:,3),temp_noent(:,7),temp_noent(:,4:6),'filled','MarkerFaceAlpha',0.8); %make no ent channels slightly transparent
            end
            %identify significant electrodes and plot:
            temp_ent=ent_mat(ent_mat(:,end)~=1,:);
            if ~isempty(temp_ent)
                scatter3(temp_ent(:,1),temp_ent(:,2),temp_ent(:,3),temp_ent(:,7),temp_ent(:,4:6),'filled');
            end
        end
        set(gca,'XTick',[], 'YTick', [],'ZTick',[]);
        colormap(flipud(autumn));
        
        if strcmp(fnames.task,'flickerfreq')
            title(i{:});
        end
        
        subplot_nber=subplot_nber+1;
    end
    
    %clean figure:
    hp4 = gca;
    hp4=hp4.Position;
    c=colorbar;
    c.TickLabels=arrayfun(@(x) num2str(x),0:2:threshold,'UniformOutput',false);
    c.TickLabels(1)={'>0'};
    c.TickLabels(end)={['\geq' num2str(threshold)]};
    c.Label.String='Fold increase in power';
    
    %add subplot labels:
    if strcmp(fnames.task,'flickerneuro')
        annotation('textbox', [0.15, 0.98, 0, 0], 'string', 'visual','FontSize',20,'FontWeight','bold');
        annotation('textbox', [0.45, 0.98, 0, 0], 'string', 'audio-visual','FontSize',20,'FontWeight','bold');
        annotation('textbox', [0.79, 0.98, 0, 0], 'string', 'audio','FontSize',20,'FontWeight','bold');
        annotation('textbox', [0, 0.82, 0, 0], 'string', '5.5Hz','FontSize',20,'FontWeight','bold');
        annotation('textbox', [0, 0.51, 0, 0], 'string', '40Hz','FontSize',20,'FontWeight','bold');
        annotation('textbox', [0, 0.22, 0, 0], 'string', '80Hz','FontSize',20,'FontWeight','bold');
    end
    
    %save figure:
    set(gcf,'Units','inches');
    screenposition = get(gcf,'Position');
    set(gcf,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print(h,[fnames.analysis_folder '/LFP/static_ent/sub-' fnames.subjectID '_ind-ent-ref' ref_method '.pdf'],'-dpdf','-fillpage');
    
    %create rotating gif:
    %create_gif(h,20,180:-1:-180,[fnames.analysis_folder '/LFP/static_ent/sub-' fnames.subjectID '_ind-ent-ref' ref_method '.gif']);
    
    close(h); %close figure
end

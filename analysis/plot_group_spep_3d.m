%Plots SPEP modulation on 3D model of subject's brain.
%2024/02/26

function plot_group_spep_3d(fnames,subjectIDs,ref_method)
    %get mni template data:
    mni=fetch_mni_anat(fnames.root_dir);
    
    %fetch anat data:
    subject=fetch_subject_data(fnames.root_dir,{fnames.subjectID},'anat');
    
    subject_nber=1;
    for sub=subjectIDs %for each subject
        %get SPEP significance, amplitude and timing values:
        subject(subject_nber).significance=readtable([fnames.root_dir '/stg-analyses/task-spep/sub-' sub{:} '/ses-01/LFP/spep_erp/ERP_pvalue_table_ref' ref_method '.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
        subject(subject_nber).amplitude=readtable([fnames.root_dir '/stg-analyses/task-spep/sub-' sub{:} '/ses-01/LFP/spep_erp/ERP_amp_table_ref' ref_method '.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
        subject(subject_nber).timing=readtable([fnames.root_dir '/stg-analyses/task-spep/sub-' sub{:} '/ses-01/LFP/spep_erp/ERP_timing_table_ref' ref_method '.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);

        subject_nber=subject_nber+1;
    end
    
    %plot SPEP response for all patients on MNI brain (1 brain for each condition):
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    %threshold=10;
    tiledlayout(1,3,'TileSpacing','compact');
    subplot_nber=1;
    for i={'V','AV','A'} %for each condition
        nexttile(subplot_nber);
        %organize all xyz electrode coordinates, spep amplitude, spep significance and timing into a n x 5 matrix:
        spep_mat=[];
        for sub_nber=1:length(subject) %for each subject
            for j=1:size(subject(sub_nber).significance,1) %for each electrode
                if subject(sub_nber).significance{j,i}<0.05 %significant electrode are larger and in color
                    dot_color=subject(sub_nber).amplitude{j,i}; %color represents amplitude of max peak
                    dot_size=10+10*subject(sub_nber).timing{j,i}; %size represents timing of max peak
                else %non-significant electrodes are smaller and grey
                    dot_color=-1;
                    dot_size=1;
                end
                temp_contact=strsplit(subject(sub_nber).significance.Properties.RowNames{j},'-');
                spep_mat(end+1,:)=[subject(sub_nber).anat.fs_coords_mni(strcmp(subject(sub_nber).anat.electrodes_info.labels,temp_contact{1}),:) dot_color dot_size]; %negative amplitudes are in grey
            end
        end

        %plot brain in 3D:
        ft_plot_mesh(mni.anat.pial,'facealpha',0.02,'edgecolor','none'); %plot MNI pial surface
        ft_plot_mesh(mni.anat.mesh_lh,'facealpha',0.1,'facecolor','cortex','edgecolor','none'); %plot left hippocampus
        ft_plot_mesh(mni.anat.mesh_rh,'facealpha',0.1,'facecolor','cortex','edgecolor','none'); %plot right hippocampus
        hold on;
        lighting gouraud;
        camlight;
        axis tight;
        if ~isempty(spep_mat)
            %plot non-significant spep electrodes:
            temp_nospep=spep_mat(spep_mat(:,4)==-1,:);
            if ~isempty(temp_nospep)
                scatter3(temp_nospep(:,1),temp_nospep(:,2),temp_nospep(:,3),temp_nospep(:,5),[0.5 0.5 0.5],'filled','MarkerFaceAlpha',0.8); %make no ent channels slightly transparent
            end
            %plot significant spep electrodes:
            temp_spep=spep_mat(spep_mat(:,4)~=-1,:);
            if ~isempty(temp_spep)
                scatter3(temp_spep(:,1),temp_spep(:,2),temp_spep(:,3),temp_spep(:,5),temp_spep(:,4),'filled');
            end
        end
        set(gca,'XTick',[], 'YTick', [],'ZTick',[]);
        colormap(flipud(autumn)); %apply color scheme
        
        subplot_nber=subplot_nber+1;
    end
    
    %clean figure:
    c=colorbar;
    c.Label.String='SPEP maximum absolute peak amplitude (\muV)';
    no_SPEP_legend=[c.Position(1) c.Position(2)-0.02 c.Position(3) c.Position(3)];
    annotation('rectangle',no_SPEP_legend,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5);
    no_SPEP_legend(1)=no_SPEP_legend(1)+no_SPEP_legend(3);
    no_SPEP_legend(3)=no_SPEP_legend(3)*7;
    annotation('textbox',no_SPEP_legend,'String','No SPEP','EdgeColor','none','FontSize',9); %SHOULD WE CHANGE THIS? TECHNICALLY GREY MEANS NEGATIVE SIGNIFICANT MODULATION, OR NON-SIGNIFICANT RESULTS
    han=axes(h,'visible','off');
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    annotation('textbox', [0.18, 0.98, 0, 0], 'string', 'visual','FontSize',20,'FontWeight','bold');
    annotation('textbox', [0.45, 0.98, 0, 0], 'string', 'audio-visual','FontSize',20,'FontWeight','bold');
    annotation('textbox', [0.75, 0.98, 0, 0], 'string', 'audio','FontSize',20,'FontWeight','bold');
    
    %save figure:
    print(h,[fnames.root_dir '/stg-analyses/task-spep/group-spep-ref-' ref_method '-' num2str(length(subjectIDs)) 'subjects-' num2str(size(spep_mat,1)) 'contacts.pdf'],'-dpdf','-fillpage');
    %create_gif(h,20,180:-1:-180,[fnames.root_dir '/stg-analyses/task-spep/group-ent-ref-' ref_method '-' num2str(length(subjectIDs)) 'subjects-' num2str(size(spep_mat,1)) 'contacts.gif']);
    close(h);
end

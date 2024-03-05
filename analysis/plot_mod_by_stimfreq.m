%Plots PSD response to flicker stimulation, for all flicker frequencies tested.
%2024/02/25

function plot_mod_by_stimfreq(fnames,ref_method)
    %fetch data:
    mod_sigs=readtable([fnames.analysis_folder '/LFP/static_ent/LFP_pvalue_table_ref' ref_method '.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true); %get modulation significance values
    mod_amps=readtable([fnames.analysis_folder '/LFP/static_ent/LFP_zscore_table_ref' ref_method '.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true); %get modulation amplitude values
    PSD_results=importdata([fnames.preproc_LFPdata '/sub-' fnames.subjectID '_stg-analysis_task-' fnames.task '_ses-' fnames.ses '_nat-psd-ref' ref_method '.mat'],'PSD_results_ref_preproc'); %get PSD results
    subject=fetch_subject_data(define_flicker_root_dir,{fnames.subjectID},'anat'); %get anatomical data
    
    %convert electrode coordinates to freesurfer space:
    temp=subject.anat.electrodes_info.voxcoords_postopct;
    temp(:,4)=1;
    temp=temp';
    temp=subject.anat.electrodes_info.postopct.vox2ras*temp;
    temp=subject.anat.electrodes_info.postopct.ras2ras*temp;
    temp=temp';
    temp(:,4)=[];
    subject.anat.electrodes_info.ecoords_native_world=temp;
    subject.anat.pial=ft_read_headshape({[fnames.root_dir '/stg-preproc/sub-' fnames.subjectID '/anat/segmentation/sub-' fnames.subjectID '_freesurfer-output/surf/lh.pial.T1'],[fnames.root_dir '/stg-preproc/sub-' fnames.subjectID '/anat/segmentation/sub-' fnames.subjectID '_freesurfer-output/surf/rh.pial.T1']});
    subject.anat.fs_t1=ft_read_mri([fnames.root_dir '/stg-preproc/sub-' fnames.subjectID '/anat/segmentation/sub-' fnames.subjectID '_freesurfer-output/mri/brain.mgz']);
    %transform to ras space:
    transform=subject.anat.fs_t1.transform*inv(subject.anat.fs_t1.hdr.tkrvox2ras);
    temp_pos=subject.anat.pial.pos;
    temp_pos(:,4)=1;
    temp_pos=temp_pos';
    subject.anat.pial.pos=transform*temp_pos;
    subject.anat.pial.pos=subject.anat.pial.pos';
    subject.anat.pial.pos(:,4)=[];

    %define output folder:
    outputFolder=[fnames.analysis_folder,'/LFP/static_ent/ent-by-stimfreq'];
    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end

    %plot modulation plots for all electrode contacts:
    depth_electrodes=extract_clinLFP_labels(mod_amps.Properties.RowNames); %organize electrode contacts by depth electrode
    view_configs=[-180 0;-90 0;0 90];
    for i=1:length(depth_electrodes) %for each depth electrode
        figure('units','normalized','outerposition',[0 0 1 1]);
        tiledlayout(6,7,'TileSpacing','none','Padding','none');

        %plot depth electrode on 3D brain:
        plot_locations=[1 8 15];
        for plot_nber=1:3
            nexttile(plot_locations(plot_nber));
            switch plot_nber
                case 1
                    title('Front view (coronal)');
                case 2
                    title('Left view (sagittal)');
                case 3
                    title('Top view (axial)');
            end
            ft_plot_mesh(subject.anat.pial,'facealpha',0.02,'edgecolor','none');
            hold on;
            view(view_configs(plot_nber,1),view_configs(plot_nber,2));
            lighting gouraud;
            camlight;
            axis tight;
            temp_coords=subject.anat.electrodes_info.ecoords_native_world(ismember(subject.anat.electrodes_info.labels,regexprep(depth_electrodes(i).channel_names,'-.+','')),:);
            scatter3(temp_coords(:,1),temp_coords(:,2),temp_coords(:,3),5,'k','filled','MarkerFaceAlpha',1);
        end
        
        %plot PSD plot for each electrode contact and condition:
        plot_locations=[2:7 9:14 16:21 23:28 30:35 37:42];
        for j=1:length(depth_electrodes(i).channel_names) %for each electrode contact
            nexttile(plot_locations(j));
            plot_PSD('Baseline',depth_electrodes(i).channel_names{j},PSD_results,PSD_results.label,PSD_results.condition,'k',1,0,1); %plot baseline
            hold on;

            for cond=mod_amps.Properties.VariableNames %for each condition
                freq_index=str2double(regexprep(cond,'Hz-.+',''));
                freq_index=[freq_index-1 freq_index+1];

                line_color=condition_color(regexprep(cond,'.+-','')); %determine condition color

                frequencies=PSD_results.data{strcmp(PSD_results.label,depth_electrodes(i).channel_names{j}),strcmp(PSD_results.condition,cond)}{3};
                freq_index=find(frequencies>=freq_index(1) & frequencies<=freq_index(2));

                %plot PSD for that condition:
                psd_result=PSD_results.data{strcmp(PSD_results.label,depth_electrodes(i).channel_names{j}),strcmp(PSD_results.condition,cond)}{1};
                plot(PSD_results.data{strcmp(PSD_results.label,depth_electrodes(i).channel_names{j}),strcmp(PSD_results.condition,cond)}{3}(freq_index),log10(mean(psd_result(:,freq_index))),'Color',line_color);

                %plot standard error:
                x=[PSD_results.data{strcmp(PSD_results.label,depth_electrodes(i).channel_names{j}),strcmp(PSD_results.condition,cond)}{3}(freq_index) PSD_results.data{strcmp(PSD_results.label,depth_electrodes(i).channel_names{j}),strcmp(PSD_results.condition,cond)}{3}(freq_index(end:-1:1))];
                p=patch(x,log10([mean(psd_result(:,freq_index))-std(psd_result(:,freq_index))/sqrt(size(psd_result,1)) mean(psd_result(:,freq_index(end:-1:1)))+std(psd_result(:,freq_index(end:-1:1)))/sqrt(size(psd_result,1))]),line_color,'FaceAlpha',0.2,'EdgeColor','none');
                %p=patch(x,log10([mean(psd_result)-std(psd_result) mean(psd_result(:,end:-1:1))+std(psd_result(:,end:-1:1))]),color,'FaceAlpha',0.2,'EdgeColor','none');
                set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

            end
            
            %clean plot:
            xlim([0 85.5]);
            ylabel('Power (log_1_0)');
            title(depth_electrodes(i).channel_names{j});
        end
        
        %plot amplitude of modulation by frequency of stimulation:
        for j=1:length(depth_electrodes(i).channel_names)
            nexttile(plot_locations(18+j));
            plot_ent_by_freq(mod_amps,mod_sigs,depth_electrodes(i).channel_names{j},1);
        end

        %save figure:
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,...
            'PaperPosition',[0 0 screenposition(3:4)],...
            'PaperSize',[screenposition(3:4)]);
        print(gcf,[outputFolder '/depth-electrode-' depth_electrodes(i).depth_electrode_name '_ent-by-stimfreq.pdf'],'-dpdf','-fillpage');
        close;
    end
end

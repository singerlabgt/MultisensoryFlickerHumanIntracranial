%Runs PSD analysis for flicker trials, using the Chronux toolbox.
%THINK ABOUT ADDING STEP TO GET RID OF ARTIFACTS
%2024/02/25

function static_ent_PSDanalysis(fnames)
    %Define output folder and load trials data:
    outputFolder=[fnames.preprocdata_folder '/LFP/static_ent'];
    disp('Loading trials timestamp data...')
    load(fnames.trials_filename);

    %create output folder if it does not exist:
    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end
    
    for ref_method={'Laplacian'}
        %load preprocessed LFP data:
        load([fnames.preproc_LFPdata '/sub-' fnames.subjectID '_stg-preproc_task-' fnames.task '_ses-' fnames.ses '_nat-trialsdata-ref' ref_method{:} '-preproc.mat'],'data_ref_preproc');
        
        %run PSD analysis:
        disp('Running PSD analysis...');
        PSD_results_ref_preproc=run_PSD_flicker(data_ref_preproc,trials);
        PSD_results_ref_preproc.ref=ref_method;

        %plot and save summary PSD plots for all channels, organized by depth electrode:
        if ~exist([outputFolder '/entrainment-PSDs-' ref_method{:}],'dir')
            mkdir([outputFolder '/entrainment-PSDs-' ref_method{:}]);
        end
        plot_depthelectrode_PSDs(PSD_results_ref_preproc,[outputFolder '/entrainment-PSDs-' ref_method{:}]);
        
        %plot and save summary PSD plots for the occluded condition(s), for all channels, organized by depth electrode:
        if ~exist([outputFolder '/occluded-PSDs-' ref_method{:}],'dir')
            mkdir([outputFolder '/occluded-PSDs-' ref_method{:}]);
        end
        plot_depthelectrode_occluded(PSD_results_ref_preproc,[outputFolder '/occluded-PSDs-' ref_method{:}]);

        %save results and clear data:
        disp('Saving PSD results...');
        save([outputFolder '/sub-' fnames.subjectID '_stg-analysis_task-' fnames.task '_ses-' fnames.ses '_nat-psd-ref' ref_method{:} '.mat'],'PSD_results_ref_preproc','-v7.3');
        clear data_ref_preproc PSD_results_ref_preproc;
    end
end

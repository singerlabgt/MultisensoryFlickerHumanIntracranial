%Fetches analyzed data from subjects.
%NEED TO PROOFCHECK
%2024/02/26

function subject=fetch_subject_data(root_dir,subjectIDs,varargin)
    subject=struct(); %initialize structure in which will store values
    
    for sub_nber=1:length(subjectIDs) %for each subject
        subject(sub_nber).subjectID=subjectIDs{sub_nber}; %subject ID
        
        %fetch anat data if requested:
        if any(ismember(varargin,'anat'))
            subject(sub_nber).anat.fs_outputdir=[root_dir '/stg-preproc/sub-' subject(sub_nber).subjectID '/anat/segmentation/sub-' subject(sub_nber).subjectID '_freesurfer-output'];
            %subject(sub_nber).anat.pial=ft_read_headshape({[subject(sub_nber).anat.fs_outputdir '/surf/lh.pial.T1'],[subject(sub_nber).anat.fs_outputdir '/surf/rh.pial.T1']});
            subject(sub_nber).anat.electrodes_info=importdata([root_dir '/stg-preproc/sub-' subject(sub_nber).subjectID '/anat/electrode-models/sub-' subject(sub_nber).subjectID '_electrodes_info.mat'],'electrodes_info');
            subject(sub_nber).anat.fs_coords_mni=subject(sub_nber).anat.electrodes_info.ecoords_mni_world;
        end
        
        %fetch task data:
        tasks=varargin(~ismember(varargin,'anat'));
        for t=tasks %for each task
            temp=strsplit(t{:},':');
            task_name=temp{1}; %task
            set_names=temp{2}; %requested results from task
            
            % check whether there are multiple sessions for this subject and task:
            temp=getElementPaths([root_dir '/stg-preproc/sub-' subject(sub_nber).subjectID '/task-' task_name],'ses-*','folder'); %get paths to folders which name starts with 'ses-'
            temp=regexprep(temp,'.+(\\|/)',''); %get just the names of those folders (instead of whole paths)
            if length(temp)>1 && contains(t,'flickerneuro') %means that we have more than 1 session for flickerneuro- need to pick sesson with strongest modulation at 40Hz-AV
                disp(['For ' subject(sub_nber).subjectID ', detected ' num2str(length(temp)) ' sessions...']);
                for i=1:length(temp) %for each session
                    temp_pvalue=readtable([root_dir '/stg-analyses/task-' task_name '/sub-' subject(sub_nber).subjectID '/ses-' regexprep(temp{i},'ses-','') '/LFP/static_ent/LFP_pvalue_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true); %get modulation significance values
                    temp_pvalue=temp_pvalue.('40Hz-AV'); %only keep results for the 40Hz-AV condition
                    temp_zscore{i}=readtable([root_dir '/stg-analyses/task-' task_name '/sub-' subject(sub_nber).subjectID '/ses-' regexprep(temp{i},'ses-','') '/LFP/static_ent/LFP_zscore_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true); %get modulation amplitude values
                    temp_zscore{i}=temp_zscore{i}.('40Hz-AV'); %only keep results for 40Hz-AV condition
                    temp_zscore{i}=temp_zscore{i}.*(temp_pvalue<=0.05); %only keep results that are significant
                end
                [~,session_nber]=max(arrayfun(@(x) sum(x{:}>0),temp_zscore)); %find session with highest number of contacts entrained in 40Hz-AV condition
                disp(['Picking ' char(temp(session_nber)) ' because highest number of significantly modulated contacts at 40Hz-AV.']);
                session_name={char(temp(session_nber))};
            elseif length(temp)>1 && contains(t,'flickerfreq') %means that we have more than 1 session for the flickerfreq task- will need to retrieve all sessions
                session_name=arrayfun(@(x) char(x{:}),temp,'UniformOutput',false);
            else %means we only have 1 session
                session_name={char(temp)};
            end
            
            %fetch results:
            set_names=strsplit(set_names,','); %requested results from task
            for set_name=set_names %for each requested result
                switch set_name{:}
                    case 'ssep_amp' %flicker modulation results
                        if length(session_name)>1 %if more than 1 session, retrieve data from all sessions
                            %initialize cell arrays for results:
                            subject(sub_nber).([task_name '_ssep_amp_sig'])=cell(length(session_name),1);
                            subject(sub_nber).([task_name '_ssep_amp_val'])=cell(length(session_name),1);
                            %fetch modulation significance and amplitude values:
                            for s=1:length(session_name)
                                subject(sub_nber).([task_name '_ssep_amp_sig']){s}=readtable([root_dir '/stg-analyses/task-' task_name '/sub-' subject(sub_nber).subjectID '/' session_name{s} '/LFP/static_ent/LFP_pvalue_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
                                subject(sub_nber).([task_name '_ssep_amp_val']){s}=readtable([root_dir '/stg-analyses/task-' task_name '/sub-' subject(sub_nber).subjectID '/' session_name{s} '/LFP/static_ent/LFP_zscore_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
                            end
                        else
                            subject(sub_nber).([task_name '_ssep_amp_sig'])=readtable([root_dir '/stg-analyses/task-' task_name '/sub-' subject(sub_nber).subjectID '/' session_name{:} '/LFP/static_ent/LFP_pvalue_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
                            subject(sub_nber).([task_name '_ssep_amp_val'])=readtable([root_dir '/stg-analyses/task-' task_name '/sub-' subject(sub_nber).subjectID '/' session_name{:} '/LFP/static_ent/LFP_zscore_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
                        end
                    case 'ssep_plv' %flicker PLV results
                        if length(session_name)>1 %if more than 1 session, retrieve data from all sessions
                            %initialize cell arrays for results:
                            subject(sub_nber).([task_name '_ssep_plv_sig'])=cell(length(session_name),1);
                            subject(sub_nber).([task_name '_ssep_plv_amp'])=cell(length(session_name),1);
                            %fetch PLV significance and amplitude values:
                            for s=1:length(session_name)
                                subject(sub_nber).([task_name '_ssep_plv_sig']){s}=readtable([root_dir '/stg-analyses/task-' task_name '/sub-' subject(sub_nber).subjectID '/' session_name{s} '/LFP/static_ent/LFP_plv-pval_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
                                subject(sub_nber).([task_name '_ssep_plv_amp']){s}=readtable([root_dir '/stg-analyses/task-' task_name '/sub-' subject(sub_nber).subjectID '/' session_name{s} '/LFP/static_ent/LFP_plv-amp_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
                            end
                        else
                            subject(sub_nber).([task_name '_ssep_plv_sig'])=readtable([root_dir '/stg-analyses/task-' task_name '/sub-' subject(sub_nber).subjectID '/' session_name{:} '/LFP/static_ent/LFP_plv-pval_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
                            subject(sub_nber).([task_name '_ssep_plv_amp'])=readtable([root_dir '/stg-analyses/task-' task_name '/sub-' subject(sub_nber).subjectID '/' session_name{:} '/LFP/static_ent/LFP_plv-amp_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true);
                        end
                    case 'spep_amp' %spep results
                        subject(sub_nber).([task_name '_spep_sig'])=readtable([root_dir '/stg-analyses/task-' task_name '/sub-' subject(sub_nber).subjectID '/' session_name{:} '/LFP/spep_erp/ERP_pvalue_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true); %fetch spep significance values
                        subject(sub_nber).([task_name '_spep_amp_val'])=readtable([root_dir '/stg-analyses/task-' task_name '/sub-' subject(sub_nber).subjectID '/' session_name{:} '/LFP/spep_erp/ERP_amp_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true); %fetch spep amplitude values
                        subject(sub_nber).([task_name '_spep_timing_val'])=readtable([root_dir '/stg-analyses/task-' task_name '/sub-' subject(sub_nber).subjectID '/' session_name{:} '/LFP/spep_erp/ERP_timing_table_refLaplacian.csv'],'ReadVariableNames',true,'ReadRowNames',true,'PreserveVariableNames',true); %fetch spep timing values
                    case 'psd_data' %PSD results
                        if length(session_name)>1 %if more than 1 session, retrieve data from all sessions
                            %initialize cell array for results
                            subject(sub_nber).([task_name '_psd_data'])=cell(length(session_name),1);
                            for s=1:length(session_name)
                                subject(sub_nber).([task_name '_psd_data']){s}=importdata([root_dir '/stg-preproc/sub-' subject(sub_nber).subjectID '/task-' task_name '/' session_name{s} '/LFP/static_ent/sub-' subject(sub_nber).subjectID '_stg-analysis_task-' task_name '_' session_name{s} '_nat-psd-refLaplacian.mat'],'PSD_results_ref_preproc');
                            end
                        else
                            subject(sub_nber).([task_name '_psd_data'])=importdata([root_dir '/stg-preproc/sub-' subject(sub_nber).subjectID '/task-' task_name '/' session_name{:} '/LFP/static_ent/sub-' subject(sub_nber).subjectID '_stg-analysis_task-' task_name '_' session_name{:} '_nat-psd-refLaplacian.mat'],'PSD_results_ref_preproc');
                        end
                end
            end
        end
    end
end

%Estimates the significance and amplitude of flicker modulation for each channel and condition, and saves those estimates into tables.
%2024/02/25

function evaluate_ent_degree_relpower(fnames)
    %create output folder if it does not exist
    if ~exist([fnames.analysis_folder,'/LFP/static_ent'],'dir')
        mkdir([fnames.analysis_folder,'/LFP/static_ent']);
    end
    
    %quantify significance and amplitude of modulation:
    for ref_method={'Laplacian'}
        %load PSD results:
        load([fnames.preproc_LFPdata '/sub-' fnames.subjectID '_stg-analysis_task-' fnames.task '_ses-' fnames.ses '_nat-psd-ref' ref_method{:} '.mat'],'PSD_results_ref_preproc');
        PSD_results=PSD_results_ref_preproc; %choose which version  of PSD results you want to use

        %define conditions we are interested in looking at:
        conditions_of_interest=sort(PSD_results.condition(contains(PSD_results.condition,'Hz') & ~contains(PSD_results.condition,'occluded') & ~contains(PSD_results.condition,'min')));
        control_condition='Baseline'; %define baseline condition
        
        %calculate modulation significance and amplitude for each channel and condition:
        pvalue_table=zeros(length(PSD_results.label),length(conditions_of_interest)); %initialize significance table
        zscore_table=zeros(length(PSD_results.label),length(conditions_of_interest)); %initialize amplitude table
        for i=1:size(zscore_table,2) %for each condition
            freq_interest=str2num(regexprep(conditions_of_interest{i},'Hz.+','')); %find frequency of stimulation
            [~,temp]=min(abs(PSD_results.data{1,strcmp(PSD_results.condition,conditions_of_interest{i})}{3}-freq_interest)); %find index of PSD analysis frequency closest to frequency of interest- have to do this because sample rate of EDF file not an integer sometimes (error with Natus)
            freq_interest_index=temp;
            
            for ch=1:size(zscore_table,1) %for each channel
                stim_values=[];
                baseline_values=[];
                
                for tr=1:size(PSD_results.data{ch,strcmp(PSD_results.condition,conditions_of_interest{i})}{1},1) %for however many trials there are for given condition
                    current_stim_value=PSD_results.data{ch,strcmp(PSD_results.condition,conditions_of_interest{i})}{1}(tr,freq_interest_index); %get PSD value at frequency of stimulation, in the stimulation trial
                    current_baseline_value=PSD_results.data{ch,strcmp(PSD_results.condition,control_condition)}{1}(tr,freq_interest_index); %get PSD value at frequency of stimulation, in the baseline trial

                    %concatenate stim and baseline values:
                    stim_values=[stim_values current_stim_value];
                    baseline_values=[baseline_values current_baseline_value];
                end
                
                %calculate modulation significance and amplitude, based on trial values:
                pvalue_table(ch,i)=pval_randomshuffle([stim_values' baseline_values'],10000);
                zscore_table(ch,i)=(mean(stim_values)/mean(baseline_values))-1;
            end
        end

        %save modulation significance and amplitude values into tables:
        pvalue_table=array2table(pvalue_table,'RowNames',PSD_results.label,'VariableNames',conditions_of_interest); %make matrix into table
        writetable(pvalue_table,[fnames.analysis_folder,'/LFP/static_ent/LFP_pvalue_table_ref' ref_method{:} '.csv'],'WriteRowNames',1); %write table
        zscore_table=array2table(zscore_table,'RowNames',PSD_results.label,'VariableNames',conditions_of_interest); %make matrix into table
        writetable(zscore_table,[fnames.analysis_folder,'/LFP/static_ent/LFP_zscore_table_ref' ref_method{:} '.csv'],'WriteRowNames',1); %write table
        
        %clear variables:
        clear PSD_results_ref_preproc PSD_results;
    end
end

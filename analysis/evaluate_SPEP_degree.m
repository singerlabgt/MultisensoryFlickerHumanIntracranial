%Evaluates significance and amplitude of SPEP.
%2024/02/25

function evaluate_SPEP_degree(fnames)
    %create output folder if it does not exist:
    if ~exist([fnames.analysis_folder,'/LFP/spep_erp'],'dir')
        mkdir([fnames.analysis_folder,'/LFP/spep_erp']);
    end

    for ref_method={'Laplacian'}
        ERP_results=importdata([fnames.preproc_LFPdata '/sub-' fnames.subjectID '_stg-analysis_task-' fnames.task '_ses-' fnames.ses '_nat-spepERP-ref' ref_method{:} '.mat'],'ERP_results_ref_preproc'); %fetch ERP results
        
        conditions_of_interest=["V";"AV";"A"]; %define conditions of interest
        %control_condition='occluded_AV';
        pvalue_table=zeros(length(ERP_results{1}.label),length(conditions_of_interest)); %will store significance of SPEP response
        amp_table=zeros(length(ERP_results{1}.label),length(conditions_of_interest)); %will store amplitude of absolute max peak
        timing_table=zeros(length(ERP_results{1}.label),length(conditions_of_interest)); %will store time of absolute max peak
        time_interest_index=ERP_results{1}.time>=0 & ERP_results{1}.time<=1; %find index of time from start of pulse to 1s after start of pulse
        time_interest=ERP_results{1}.time(time_interest_index);
        for i=1:size(pvalue_table,2) %for each condition
            for ch=1:size(pvalue_table,1) %for each channel
%                 stim_values=[];
%                 baseline_values=[];
%                 
%                 if size(ERP_results{i}.trial,1)~=size(ERP_results{4}.trial,1) %check that stim and occluded conditions have same number of trials
%                     error(['Not same number of trials for ' conditions_of_interest{i} ' vs occluded.']);
%                 end
%                 
%                 for tr=1:size(ERP_results{i}.trial,1) %for however many number of trials of given condition there are
%                     current_stim_value=rms(squeeze(ERP_results{i}.trial(tr,ch,time_interest_index)));
%                     current_baseline_value=rms(squeeze(ERP_results{4}.trial(tr,ch,time_interest_index)));
%                     
%                     stim_values=[stim_values current_stim_value];
%                     baseline_values=[baseline_values current_baseline_value];
%                 end
%                 pvalue_table(ch,i)=pval_randomshuffle([stim_values' baseline_values'],10000);
                
                stim_signal=reshape(ERP_results{i}.trial(:,ch,time_interest_index),size(ERP_results{i}.trial,1),sum(time_interest_index));
                average_stim_signal=mean(stim_signal);
                
                [max_peak,time_max_peak]=max(abs(average_stim_signal));
                amp_table(ch,i)=max_peak;
                timing_table(ch,i)=time_interest(time_max_peak);
                
                %trying a different method:
                average_baseline_signal=mean(reshape(ERP_results{4}.trial(:,ch,time_interest_index),size(ERP_results{4}.trial,1),sum(time_interest_index)));
                difference_rms=rms(average_stim_signal)-rms(average_baseline_signal);
                
                x=0;
                num_shuffles=500;
                stim_and_baseline_trials=[reshape(ERP_results{i}.trial(:,ch,time_interest_index),size(ERP_results{i}.trial,1),sum(time_interest_index));reshape(ERP_results{4}.trial(:,ch,time_interest_index),size(ERP_results{4}.trial,1),sum(time_interest_index))];
                temp_different_rms=0;
                length_trials=size(stim_and_baseline_trials,1);
                length_half_trials=length_trials/2;
                for j=1:num_shuffles
                    temp_index=randperm(length_trials);
                    temp_different_rms=rms(mean(stim_and_baseline_trials(temp_index(1:length_half_trials),:)))-rms(mean(stim_and_baseline_trials(temp_index(length_half_trials+1:end),:)));
                    if temp_different_rms>difference_rms
                        x=x+1;
                    end
                end
                pvalue_table(ch,i)=x/num_shuffles;
                
            end
        end
        
        %save results:
        pvalue_table=array2table(pvalue_table,'RowNames',ERP_results{1}.label,'VariableNames',conditions_of_interest); %make matrix into table
        writetable(pvalue_table,[fnames.analysis_folder,'/LFP/spep_erp/ERP_pvalue_table_ref' ref_method{:} '.csv'],'WriteRowNames',1);
        
        amp_table=array2table(amp_table,'RowNames',ERP_results{1}.label,'VariableNames',conditions_of_interest); %make matrix into table
        writetable(amp_table,[fnames.analysis_folder,'/LFP/spep_erp/ERP_amp_table_ref' ref_method{:} '.csv'],'WriteRowNames',1);
        
        timing_table=array2table(timing_table,'RowNames',ERP_results{1}.label,'VariableNames',conditions_of_interest); %make matrix into table
        writetable(timing_table,[fnames.analysis_folder,'/LFP/spep_erp/ERP_timing_table_ref' ref_method{:} '.csv'],'WriteRowNames',1);
        
%         amp_table=(table2array(pvalue_table)<0.05).*amp_table{:,:};
%         timing_table=(table2array(pvalue_table)<0.05).*timing_table{:,:};
%         
%         amp_table=array2table(amp_table,'RowNames',ERP_results{1}.label,'VariableNames',conditions_of_interest); %make matrix into table
%         timing_table=array2table(timing_table,'RowNames',ERP_results{1}.label,'VariableNames',conditions_of_interest); %make matrix into table
        
        clear ERP_results_ref_preproc ERP_results;
    end
end

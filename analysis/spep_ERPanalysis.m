%Performs ERP analysis on spep preprocessed LFP.
%2024/02/25

function spep_ERPanalysis(fnames)
    
    %define filenames and load data:
    outputFolder=[fnames.preprocdata_folder '/LFP/spep_erp'];
    disp('Loading trials timestamp data...')
    load(fnames.trials_filename);
    
    %run ERP analysis:
    if ~exist(outputFolder,'dir')
        mkdir(outputFolder);
    end
    
    for ref_method={'Laplacian'}
        disp('Running ERP analysis...');
        data=importdata([fnames.preproc_LFPdata '/sub-' fnames.subjectID '_stg-preproc_task-' fnames.task '_ses-' fnames.ses '_nat-trialsdata-ref' ref_method{:} '-preproc.mat'],'data_ref_preproc'); %fetch preprocessed data

        conditions_to_test=["V";"AV";"A";"occluded_AV"];
        %make cell array of ERP results for 10s trials:
        %neeed to use FieldTrip's 'nearest' function here, otherwise error.
        ERP_results_ref_preproc=cell(1,length(conditions_to_test));
        for condition=conditions_to_test' %for each condition
            cfg=[];
            cfg.trials=find(data.trialinfo==find(strcmp(trials.condition_code,condition)));
            cfg.keeptrials='yes';
            ERP_results_ref_preproc{strcmp(conditions_to_test,condition)}=ft_timelockanalysis(cfg,data);
        end
        
        ERP_results_ref_preproc{end+1}=['ref_method: ' ref_method{:}];
        ERP_results_ref_preproc{end+1}=conditions_to_test;
  
        %plot results:
        disp('Plotting results...');
        if ~exist([outputFolder '/spep-ERPs-' ref_method{:}],'dir')
            mkdir([outputFolder '/spep-ERPs-' ref_method{:}]);
        end
        conditions_of_interest=conditions_to_test;
        time_to_plot=[-0.25 1]; %in s
        [~,index1]=min(abs(ERP_results_ref_preproc{1,1}.time-time_to_plot(1)));
        [~,index2]=min(abs(ERP_results_ref_preproc{1,1}.time-time_to_plot(2)));
        line_width=2;
        clinLFP_contacts=extract_clinLFP_labels(ERP_results_ref_preproc{1,1}.label);
        fig=[];
        for depth_electrode=1:length(clinLFP_contacts) %for each depth electrode
            disp(['Plotting depth electrode ' num2str(depth_electrode) ' out of ' num2str(length(clinLFP_contacts)) ' depth electrodes.']);
            fig(depth_electrode)=figure('units','normalized','outerposition',[0 0 1 1]);
            tiledlayout(4,5,'Padding','none','TileSpacing','none');
            for ch=clinLFP_contacts(depth_electrode).channel_names' %for each electrode contact
                nexttile;
                hold on;
                for current_condition=conditions_of_interest' %for each condition
                    line_color=condition_color(current_condition); %determine condition color
                    
                    erp_result=reshape(ERP_results_ref_preproc{strcmp(conditions_to_test,current_condition)}.trial(:,strcmp(ERP_results_ref_preproc{strcmp(conditions_to_test,current_condition)}.label,ch),:),size(ERP_results_ref_preproc{strcmp(conditions_to_test,current_condition)}.trial,1),size(ERP_results_ref_preproc{strcmp(conditions_to_test,current_condition)}.trial,3));
                    
                    %plot ERP:
                    plot(ERP_results_ref_preproc{strcmp(conditions_to_test,current_condition)}.time,mean(erp_result),'Color',line_color,'LineWidth',line_width);
                    
                    %plot standard error:
                    x=[ERP_results_ref_preproc{strcmp(conditions_to_test,current_condition)}.time ERP_results_ref_preproc{strcmp(conditions_to_test,current_condition)}.time(end:-1:1)];
                    p=patch(x,[mean(erp_result)-std(erp_result)/sqrt(size(erp_result,1)) mean(erp_result(:,end:-1:1))+std(erp_result(:,end:-1:1))/sqrt(size(erp_result,1))],line_color,'FaceAlpha',0.2,'EdgeColor','none');
                    %p=patch(x,log10([mean(psd_result)-std(psd_result) mean(psd_result(:,end:-1:1))+std(psd_result(:,end:-1:1))]),color,'FaceAlpha',0.2,'EdgeColor','none');
                    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                
                %clean plot:
                temp1=get(gca,'ylim');
                x = [0; 0.0125; 0.0125; 0];
                y = [temp1(1); temp1(1); temp1(2); temp1(2)];
                p=patch(x,y,[0.5 0.5 0.5],'EdgeColor','none');
                set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                xline(0,'-','','LabelHorizontalAlignment','left','FontSize',25);
                hold on;
                xlim([ERP_results_ref_preproc{1,1}.time(index1) ERP_results_ref_preproc{1,1}.time(index2)]);
               
                %xlabel('Time (ms)','FontSize',28);
                %ylabel('Electric Potential (uV)','FontSize',28);
                %set(get(gca,'YLabel'),'Rotation',90);
                plot([ERP_results_ref_preproc{strcmp(conditions_to_test,current_condition)}.time(1), ERP_results_ref_preproc{strcmp(conditions_to_test,current_condition)}.time(end)], [0 0], 'k--'); % add horizontal line
                %title(['Channel ' char(ch) newline 'anat: ' clinLFP_contacts(depth_electrode).channel_names{ismember(clinLFP_contacts(depth_electrode).channel_names(:,1),ch),2}]);
                temp=strsplit(ch{:},'-');
                title(ch);
                %title(['Channel ' char(ch) newline 'anat: ' char(electrodes_info.anatlabels{strcmp(electrodes_info.labels,temp{1}),'Neuro'})],'FontSize',15);
                %legend(conditions_of_interest,'Interpreter','none','FontSize',20);
                %legend boxoff; 
                set(gca,'children',flipud(get(gca,'children')));
            end
            
            %save figure:
            set(gcf,'Units','inches');
            screenposition = get(gcf,'Position');
            set(gcf,...
                'PaperPosition',[0 0 screenposition(3:4)],...
                'PaperSize',[screenposition(3:4)]);
            print(fig(depth_electrode),[outputFolder '/spep-ERPs-' ref_method{:} '/depth-electrode-' clinLFP_contacts(depth_electrode).depth_electrode_name '_spep-ERP.pdf'],'-dpdf','-fillpage');
        end
        
        close all;
%         xlim([-0.05 0.3])
%         ax=gca;
%         temp=ax.XTickLabel;
%         ax.XTickLabel=cellfun(@(x) num2str(str2num(x)*1000),temp,'UniformOutput',false);
        
        %save results:
        disp('Saving ERP results...');
        save([outputFolder '/sub-' fnames.subjectID '_stg-analysis_task-' fnames.task '_ses-' fnames.ses '_nat-spepERP-ref' ref_method{:} '.mat'],'ERP_results_ref_preproc','-v7.3');
    end
end

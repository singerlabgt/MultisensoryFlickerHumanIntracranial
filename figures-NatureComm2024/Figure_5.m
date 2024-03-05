%
%2024/02/26

%define directories
root_dir=define_flicker_root_dir;
source_data=[root_dir '/stg-analyses/NatureComm2024-figures/Source Data']; %where source data for figures is stored
figures_dir=[root_dir '/stg-analyses/NatureComm2024-figures'];

%fetch data:
all_subjects=fetch_flicker_subjectIDs(root_dir,'flickerfreq');
subject=fetch_subject_data(root_dir,all_subjects','flickerfreq:ssep_amp,ssep_plv,psd_data','anat');

%gather all results in summary tables (using both fold-change values and plv values):
conditions=order_flicker_conditions(subject(1).flickerfreq_ssep_amp_sig.Properties.VariableNames); %order conditions
conditions=regexprep(conditions,'-.+',''); %remove information about modality
mod_summary=struct('amp_tbl',[],'plv_tbl',[]);
for tbl={'amp_tbl','plv_tbl'}
    %define fields to fetch from in subject structure:
    if strcmp(tbl,'amp_tbl')
        sig_field='flickerfreq_ssep_amp_sig';
        val_field='flickerfreq_ssep_amp_val';
    elseif strcmp(tbl,'plv_tbl')
        sig_field='flickerfreq_ssep_plv_sig';
        val_field='flickerfreq_ssep_plv_amp';
    end
    
    for sub_nber=1:length(subject)
        if iscell(subject(sub_nber).(sig_field)) %means there's multiple sessions
            for s=1:length(subject(sub_nber).(sig_field))
                temp=table2array(subject(sub_nber).(sig_field){s});
                temp(temp>0.05)=NaN;
                temp(temp<0.05)=1;
                temp=table2array(subject(sub_nber).(val_field){s}).*temp;
                temp=array2table(temp,'VariableNames',regexprep(subject(sub_nber).(sig_field){s}.Properties.VariableNames,'-.+',''),'RowNames',strcat(subject(sub_nber).subjectID,';',subject(sub_nber).(sig_field){s}.Properties.RowNames,';',regexprep(subject(sub_nber).(sig_field){s}.Properties.VariableNames{1},'.+Hz-',''),';',num2str(s)));
                temp=temp(:,conditions);
                mod_summary.(tbl{:})=[mod_summary.(tbl{:});temp];
            end
        else
            temp=table2array(subject(sub_nber).(sig_field));
            temp(temp>0.05)=NaN;
            temp(temp<0.05)=1;
            temp=table2array(subject(sub_nber).(val_field)).*temp;
            temp=array2table(temp,'VariableNames',regexprep(subject(sub_nber).(sig_field).Properties.VariableNames,'-.+',''),'RowNames',strcat(subject(sub_nber).subjectID,';',subject(sub_nber).(sig_field).Properties.RowNames,';',regexprep(subject(sub_nber).(sig_field).Properties.VariableNames{1},'.+Hz-','')));
            temp=temp(:,conditions);
            mod_summary.(tbl{:})=[mod_summary.(tbl{:});temp];
        end
    end
end

%set threshold to pick channels of interest, based on number of significant
%stimulation frequencies (i.e. channels we are confident exhibit
%resonance):
nber_sigfreq_thresh=floor(size(subject(1).flickerfreq_ssep_amp_sig,2)/4);

%plot PSDs for channels that showed modulation to more than
%nber_sigfreq_thresh stim frequencies (to check whether threshold we picked
%above is adequate):
temp=mod_summary.amp_tbl;
temp=temp(sum(~isnan(table2array(temp)'))>nber_sigfreq_thresh,:);
for i=1:size(temp,1)
    if mod(i,5*8)==1
        figure('Units','normalized','Position',[0 0 1 1]);
        tiledlayout(5,8,'Padding','none','TileSpacing','none');
    end
    nexttile;
    temp_info=strsplit(temp.Properties.RowNames{i},';');
    sub_nber=find(strcmp({subject.subjectID},temp_info{1}));
    if length(temp_info)==4
        psd_data=subject(sub_nber).flickerfreq_psd_data{str2num(temp_info{4})};
    else
        psd_data=subject(sub_nber).flickerfreq_psd_data;
    end
    
    plot_PSD('Baseline',temp_info{2},psd_data,psd_data.label,psd_data.condition,'k',1,0,1); %plot baseline PSD
    hold on;
    for cond=temp.Properties.VariableNames
        freq_index=str2double(regexprep(cond,'Hz',''));
        freq_index=[freq_index-1 freq_index+1];
        
        line_color=condition_color(regexprep(temp_info{3},'.+-',''));
        
        frequencies=psd_data.data{strcmp(psd_data.label,temp_info{2}),strcmp(psd_data.condition,[cond{:} '-' temp_info{3}])}{3};
        freq_index=find(frequencies>=freq_index(1) & frequencies<=freq_index(2));
        
        psd_result=psd_data.data{strcmp(psd_data.label,temp_info{2}),strcmp(psd_data.condition,[cond{:} '-' temp_info{3}])}{1};
        plot(psd_data.data{strcmp(psd_data.label,temp_info{2}),strcmp(psd_data.condition,[cond{:} '-' temp_info{3}])}{3}(freq_index),log10(mean(psd_result(:,freq_index))),'Color',line_color);
        
        x=[psd_data.data{strcmp(psd_data.label,temp_info{2}),strcmp(psd_data.condition,[cond{:} '-' temp_info{3}])}{3}(freq_index) psd_data.data{strcmp(psd_data.label,temp_info{2}),strcmp(psd_data.condition,[cond{:} '-' temp_info{3}])}{3}(freq_index(end:-1:1))];
        p=patch(x,log10([mean(psd_result(:,freq_index))-std(psd_result(:,freq_index))/sqrt(size(psd_result,1)) mean(psd_result(:,freq_index(end:-1:1)))+std(psd_result(:,freq_index(end:-1:1)))/sqrt(size(psd_result,1))]),line_color,'FaceAlpha',0.2,'EdgeColor','none');
        %p=patch(x,log10([mean(psd_result)-std(psd_result) mean(psd_result(:,end:-1:1))+std(psd_result(:,end:-1:1))]),color,'FaceAlpha',0.2,'EdgeColor','none');
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
    end
    
    xlim([0 85.5]);
    ylabel('Power (log_1_0)');
    
end

close all;

%diplays a few numbers for text (CHECK CODE):
disp(['Percent contacts showing significant fold change to at least 1 stim frequency: ' num2str(round(100*sum(any(~isnan(table2array(mod_summary.amp_tbl))'))/size(mod_summary.amp_tbl,1),1)) '%']);
disp(['Percent contacts showing significant fold change to more than ' num2str(nber_sigfreq_thresh) ' stim frequencies: ' num2str(round(100*sum(sum(~isnan(table2array(mod_summary.amp_tbl))')>nber_sigfreq_thresh)/size(mod_summary.amp_tbl,1),1)) '%']);
temp=mod_summary.amp_tbl;
temp=temp(sum(~isnan(table2array(temp)'))>nber_sigfreq_thresh,:);
temp=table2array(temp);
disp(['Median fold-change: ' num2str(round(median(temp(:),'omitnan'),1)), ', range: ' num2str(round(min(temp(:)),5)) ' to ' num2str(round(max(temp(:)),1)) ', 25th percentile: ' num2str(round(prctile(temp(:),25),1)) ', 75th percentile: ' num2str(round(prctile(temp(:),75),1))])
[~,temp_index]=max(temp');
disp(['Percent contacts with preference to lowest stim freq, out of contacts showing significant fold change to at least ' num2str(nber_sigfreq_thresh) ' stim frequencies: ' num2str(round(100*sum(temp_index==1)/size(temp,1),1)) '%']);
disp(['Percent contacts showing significant fold change to less than or equal to ' num2str(nber_sigfreq_thresh) ' stim frequencies: ' num2str(round(100*sum(sum(~isnan(table2array(mod_summary.amp_tbl))')<=nber_sigfreq_thresh & ~all(isnan(table2array(mod_summary.amp_tbl))'))/size(mod_summary.amp_tbl,1),1)) '%']);
temp=mod_summary.amp_tbl;
temp=temp(sum(~isnan(table2array(temp)'))<=nber_sigfreq_thresh & ~all(isnan(table2array(temp))'),:);
temp=table2array(temp);
disp(['(responding to ' num2str(nber_sigfreq_thresh) ' or less stim frequencies) Median fold-change: ' num2str(round(median(temp(:),'omitnan'),1)), ', range: ' num2str(round(min(temp(:)),1)) ' to ' num2str(round(max(temp(:)),1)) ', 25th percentile: ' num2str(round(prctile(temp(:),25),1)) ', 75th percentile: ' num2str(round(prctile(temp(:),75),1))])

disp(['Percent contacts showing significant plv to at least 1 stim frequency: ' num2str(round(100*sum(any(~isnan(table2array(mod_summary.plv_tbl))'))/size(mod_summary.plv_tbl,1),1)) '%']);
disp(['Percent contacts showing significant plv to more than ' num2str(nber_sigfreq_thresh) ' stim frequencies: ' num2str(round(100*sum(sum(~isnan(table2array(mod_summary.plv_tbl))')>nber_sigfreq_thresh)/size(mod_summary.plv_tbl,1),1)) '%']);
temp=mod_summary.plv_tbl;
temp=temp(sum(~isnan(table2array(temp)'))>nber_sigfreq_thresh,:);
temp=table2array(temp);
disp(['Median PLV: ' num2str(round(median(temp(:),'omitnan'),3)), ', range: ' num2str(round(min(temp(:)),3)) ' to ' num2str(round(max(temp(:)),3)) ', 25th percentile: ' num2str(round(prctile(temp(:),25),3)) ', 75th percentile: ' num2str(round(prctile(temp(:),75),3))])
[~,temp_index]=max(temp');
disp(['Percent contacts with preference to lowest stim freq, out of contacts showing significant plv to at least ' num2str(nber_sigfreq_thresh) ' stim frequencies: ' num2str(round(100*sum(temp_index==1)/size(temp,1),1)) '%']);
disp(['Percent contacts showing significant PLV to less than or equal to ' num2str(nber_sigfreq_thresh) ' stim frequencies: ' num2str(round(100*sum(sum(~isnan(table2array(mod_summary.plv_tbl))')<=nber_sigfreq_thresh & ~all(isnan(table2array(mod_summary.plv_tbl))'))/size(mod_summary.plv_tbl,1),1)) '%']);
temp=mod_summary.plv_tbl;
temp=temp(sum(~isnan(table2array(temp)'))<=nber_sigfreq_thresh & ~all(isnan(table2array(temp))'),:);
temp=table2array(temp);
disp(['(responding to ' num2str(nber_sigfreq_thresh) ' or less stim frequencies) Median PLV: ' num2str(round(median(temp(:),'omitnan'),3)), ', range: ' num2str(round(min(temp(:)),3)) ' to ' num2str(round(max(temp(:)),3)) ', 25th percentile: ' num2str(round(prctile(temp(:),25),3)) ', 75th percentile: ' num2str(round(prctile(temp(:),75),3))])


%% plot heatmaps of modulated channels and their modulation values:
[~,repo]=define_flicker_root_dir;
path(path,genpath(repo.fieldtrip)); %code below uses the MATLAB normalize function, as opposed to the FieldTrip normalize function; hence, puting the FieldTrip repo to end of path
iteration=0;
for tbl={'amp_tbl','plv_tbl'}
    iteration=iteration+1;
    %format data:
    modulation_matrix=table2array(mod_summary.(tbl{:}));
    modulation_matrix_labels=mod_summary.(tbl{:}).Properties.RowNames;
    
    %normalize value by channel:
    for i=1:size(modulation_matrix,1)
        modulation_matrix(i,:)=normalize(modulation_matrix(i,:),'range'); 
    end
    
    %subdivide modulation matrix into 3 groups:
    modulation_matrix_nosig=modulation_matrix(all(isnan(modulation_matrix')),:); %channels that show no response to any stim frequency
    modulation_matrix_nosig_labels=modulation_matrix_labels(all(isnan(modulation_matrix')));
    modulation_matrix_lesshalfsig=modulation_matrix(sum(~isnan(modulation_matrix'))<=nber_sigfreq_thresh & sum(~isnan(modulation_matrix'))~=0,:); %channels that show response to 13 stim frequencies or less
    modulation_matrix_lesshalfsig_labels=modulation_matrix_labels(sum(~isnan(modulation_matrix'))<=nber_sigfreq_thresh & sum(~isnan(modulation_matrix'))~=0);
    modulation_matrix_morehalfsig=modulation_matrix(sum(~isnan(modulation_matrix'))>nber_sigfreq_thresh,:); %channels that show response to 14 stim frequencies or more
    modulation_matrix_morehalfsig_labels=modulation_matrix_labels(sum(~isnan(modulation_matrix'))>nber_sigfreq_thresh);

    disp(['Percent of contacts showing response to more than ' num2str(nber_sigfreq_thresh) ' of stim frequencies: ' num2str(round(size(modulation_matrix_morehalfsig)/size(modulation_matrix)*100,2))]);
    [~,temp_index]=max(modulation_matrix_morehalfsig');
    disp(['Within, percent of contacts showing highest response at lowest stim freq: ' num2str(round(sum(temp_index==1)/length(temp_index)*100,2))]);

    % show comprehensive heatmap of types of mod-by-freq responses:
    [~,index]=max(modulation_matrix_morehalfsig');
    [~,I2]=sort(index);
    modulation_matrix_morehalfsig=modulation_matrix_morehalfsig(I2,:);
    modulation_matrix_morehalfsig_labels=modulation_matrix_morehalfsig_labels(I2);

    [~,index]=max(modulation_matrix_lesshalfsig');
    [~,I2]=sort(index);
    modulation_matrix_lesshalfsig=modulation_matrix_lesshalfsig(I2,:);
    modulation_matrix_lesshalfsig_labels=modulation_matrix_lesshalfsig_labels(I2);

    figure_making('width',3.5,'height',3);
    h=heatmap([modulation_matrix_morehalfsig; modulation_matrix_lesshalfsig],'Colormap',flipud(autumn),'GridVisible','off','MissingDataColor',[0.8 0.8 0.8],'MissingDataLabel','Not significant');
    h.FontSize=6;
    ax = gca;
    ax.XDisplayLabels = nan(size(ax.XDisplayData));
    ax.YDisplayLabels = nan(size(ax.YDisplayData));
    stim_freqs=num2cell(arrayfun(@(x) str2num(x{:}),regexprep(mod_summary.(tbl{:}).Properties.VariableNames,'Hz','')));
    x_labels=repmat({'  .  '},1,length(stim_freqs));
    x_labels(1:floor(26/5):26)=stim_freqs(1:floor(26/5):26);
    ax.XDisplayLabels=x_labels;
    xlabel('Increasing stim frequency');
    ylabel('Channels');
    
    %annotate where separation between 2 groups of channels is:
    temp=gca;
    annotation('line',[temp.InnerPosition(1)-temp.InnerPosition(1)/2 temp.InnerPosition(:,1)],repmat(temp.InnerPosition(2)+temp.InnerPosition(4)*(size(modulation_matrix_lesshalfsig,1)/(size(modulation_matrix_lesshalfsig,1)+size(modulation_matrix_morehalfsig,1))),1,2),'Color',[0.5 0.5 0.5 0.5]);

    %(size(modulation_matrix_lesshalfsig,1)+size(modulation_matrix_morehalfsig,1))/(size(modulation_matrix_nosig,1)+size(modulation_matrix_lesshalfsig,1)+size(modulation_matrix_morehalfsig,1))
    
    
    figure_making('operation','save','filename',[figures_dir '/Figure_5C_' num2str(iteration) '.pdf']);
end


%% detect endogenous oscillations for all channels and sessions:

pyenv('Version','C:/Users/lblanpa/Anaconda3/envs/fooof_env/python'); %this environmentuses python version 3.6.13 (fooof does not work with version 3.9)

plt_figure=0; %whether to plot all detected endogenous oscillations (to double-check performance of algorithm)
plt_nber=0;
for sub_nber=1:length(subject)
        
    %keep only modulation amplitudes that are significant:
    if iscell(subject(sub_nber).flickerfreq_ssep_amp_sig)
        for s=1:length(subject(sub_nber).flickerfreq_ssep_amp_sig)
            subject(sub_nber).mod_peaks{s}=subject(sub_nber).flickerfreq_ssep_amp_sig{s}; %initialize table that keeps track of significant modulation amplitudes
            temp=table2array(subject(sub_nber).mod_peaks{s});
            temp(temp>=0.05)=NaN; %make non-significant values NaN
            temp(~isnan(temp))=1; %make significant value 1
            subject(sub_nber).mod_peaks{s}{:,:}=temp.*table2array(subject(sub_nber).flickerfreq_ssep_amp_val{s}); %keep significant modulation values
            
            %retrieve endogenous freq:
            for ch=1:length(subject(sub_nber).flickerfreq_psd_data{s}.label) %for each channel
                
                %find endogenous freqs at baseline:
                
                %use FOOOF to find to 5 endogenous peaks (NEED TO FINE-TUNE):
                freqs = subject(sub_nber).flickerfreq_psd_data{s}.data{ch,strcmp(subject(sub_nber).flickerfreq_psd_data{s}.condition,'Baseline')}{3};
                psd = mean(subject(sub_nber).flickerfreq_psd_data{s}.data{ch,strcmp(subject(sub_nber).flickerfreq_psd_data{s}.condition,'Baseline')}{1});
                psd(freqs>=59 & freqs<=61)=interp1([freqs(abs(freqs-59)==min(abs(freqs-59))) freqs(abs(freqs-61)==min(abs(freqs-61)))],[psd(abs(freqs-59)==min(abs(freqs-59))) psd(abs(freqs-61)==min(abs(freqs-61)))],freqs(freqs>=59 & freqs<=61)); %draw linear interpolation at 60Hz to remove ground noise (so not picked as endogenous oscillation)
                
                settings = struct(); %initialize settings structure
                settings.max_n_peaks=5;
                settings.peak_width_limits=[2 10];
                settings.min_peak_height=0.6;
                f_range = [2, 100];
                subject(sub_nber).fooof_results{s}{ch} = fooof(freqs, psd, f_range, settings,true);
                
                %find relative change in power compared to ap_fit, and store freq, relative change in power, and range of peak:
                subject(sub_nber).endog_peaks{s}{ch}=subject(sub_nber).fooof_results{s}{ch}.peak_params;
                for i=1:size(subject(sub_nber).endog_peaks{s}{ch},1) %for each identified peak
                    [~,temp_freq_index]=min(abs(subject(sub_nber).fooof_results{s}{ch}.freqs-subject(sub_nber).endog_peaks{s}{ch}(i,1))); %find closest power spectrum frequency
                    subject(sub_nber).endog_peaks{s}{ch}(i,2)=((10^subject(sub_nber).fooof_results{s}{ch}.fooofed_spectrum(temp_freq_index))/(10^subject(sub_nber).fooof_results{s}{ch}.ap_fit(temp_freq_index)))-1; %calculate relative change in power (of model compared to ap_fit, at center frequency of endogenous oscillation) for that peak
                    if subject(sub_nber).endog_peaks{s}{ch}(i,2)<0
                        error('One of the endogenous oscillations fold change in power is negative');
                    end
                end
                
                if plt_figure && ~isempty(subject(sub_nber).endog_peaks{s}{ch})
                    plt_nber=plt_nber+1;
                    if mod(plt_nber,6*8)==1
                        figure('Units','normalized','Position',[0 0 1 1]);
                        tiledlayout(6,8,'Padding','none','TileSpacing','none');
                    end
                    nexttile;
                    plot(subject(sub_nber).fooof_results{s}{ch}.freqs,subject(sub_nber).fooof_results{s}{ch}.power_spectrum,'k');
                    hold on;
                    plot(subject(sub_nber).fooof_results{s}{ch}.freqs,subject(sub_nber).fooof_results{s}{ch}.ap_fit,'Color',[0.5 0.5 0.5]);
                    plot(subject(sub_nber).fooof_results{s}{ch}.freqs,subject(sub_nber).fooof_results{s}{ch}.fooofed_spectrum,'r');
                    for i=1:size(subject(sub_nber).endog_peaks{s}{ch},1)
                        xline(subject(sub_nber).endog_peaks{s}{ch}(i,1),'--');
                    end
                end
                
            end
        end
    else
        subject(sub_nber).mod_peaks=subject(sub_nber).flickerfreq_ssep_amp_sig; %initialize table that keeps track of significant modulation amplitudes
        temp=table2array(subject(sub_nber).mod_peaks);
        temp(temp>=0.05)=NaN; %make non-significant values NaN
        temp(~isnan(temp))=1; %make significant value 1
        subject(sub_nber).mod_peaks{:,:}=temp.*table2array(subject(sub_nber).flickerfreq_ssep_amp_val); %keep significant modulation values
        
        %retrieve endogenous freq:
        for ch=1:length(subject(sub_nber).flickerfreq_psd_data.label) %for each channel
            
            %find endogenous freqs at baseline:
            
            %use FOOOF to find to 5 endogenous peaks (NEED TO FINE-TUNE):
            freqs = subject(sub_nber).flickerfreq_psd_data.data{ch,strcmp(subject(sub_nber).flickerfreq_psd_data.condition,'Baseline')}{3};
            psd = mean(subject(sub_nber).flickerfreq_psd_data.data{ch,strcmp(subject(sub_nber).flickerfreq_psd_data.condition,'Baseline')}{1});
            psd(freqs>=59 & freqs<=61)=interp1([freqs(abs(freqs-59)==min(abs(freqs-59))) freqs(abs(freqs-61)==min(abs(freqs-61)))],[psd(abs(freqs-59)==min(abs(freqs-59))) psd(abs(freqs-61)==min(abs(freqs-61)))],freqs(freqs>=59 & freqs<=61)); %draw linear interpolation at 60Hz to remove ground noise (so not picked as endogenous oscillation)
            
            settings = struct(); %initialize settings structure
            settings.max_n_peaks=5;
            settings.peak_width_limits=[2 10];
            settings.min_peak_height=0.6;
            f_range = [2, 100];
            subject(sub_nber).fooof_results{ch} = fooof(freqs, psd, f_range, settings,true);
            %         fooof_plot(subject(sub_nber).fooof_results{ch});
            %         for i=1:size(subject(sub_nber).fooof_results{ch}.peak_params,1)
            %             xline(subject(sub_nber).fooof_results{ch}.peak_params(i,1),':','LineWidth',1);
            %         end
            
            %find relative change in power compared to ap_fit, and store freq, relative change in power, and range of peak:
            subject(sub_nber).endog_peaks{ch}=subject(sub_nber).fooof_results{ch}.peak_params;
            for i=1:size(subject(sub_nber).endog_peaks{ch},1) %for each identified peak
                [~,temp_freq_index]=min(abs(subject(sub_nber).fooof_results{ch}.freqs-subject(sub_nber).endog_peaks{ch}(i,1))); %find closest power spectrum frequency
                subject(sub_nber).endog_peaks{ch}(i,2)=((10^subject(sub_nber).fooof_results{ch}.fooofed_spectrum(temp_freq_index))/(10^subject(sub_nber).fooof_results{ch}.ap_fit(temp_freq_index)))-1;
                if subject(sub_nber).endog_peaks{ch}(i,2)<0
                    error('One of the endogenous oscillations fold change in power is negative');
                end
            end

            if plt_figure && ~isempty(subject(sub_nber).endog_peaks{ch})
                plt_nber=plt_nber+1;
                if mod(plt_nber,6*8)==1
                    figure('Units','normalized','Position',[0 0 1 1]);
                    tiledlayout(6,8,'Padding','none','TileSpacing','none');
                end
                nexttile;
                plot(subject(sub_nber).fooof_results{ch}.freqs,subject(sub_nber).fooof_results{ch}.power_spectrum,'k');
                hold on;
                plot(subject(sub_nber).fooof_results{ch}.freqs,subject(sub_nber).fooof_results{ch}.ap_fit,'Color',[0.5 0.5 0.5]);
                plot(subject(sub_nber).fooof_results{ch}.freqs,subject(sub_nber).fooof_results{ch}.fooofed_spectrum,'r');
                for i=1:size(subject(sub_nber).endog_peaks{ch},1)
                    xline(subject(sub_nber).endog_peaks{ch}(i,1),'--');
                end
            end
            %         %find peak stim freq:
            %         if any(temp_mod_amp{ch,:})
            %             [~,stim_freq]=max(temp_mod_amp{ch,:});
            %             stim_freq=temp_mod_amp.Properties.VariableNames{stim_freq};
            %             stim_freq=str2double(regexprep(stim_freq,'Hz-.+',''));
            %             subject(sub_nber).entrainment_corr{ch,2}=stim_freq;
            %         else
            %             subject(sub_nber).entrainment_corr{ch,2}=0;
            %         end
            
        end
        
        %     %convert cell array to table:
        %     subject(sub_nber).entrainment_corr=cell2table(subject(sub_nber).entrainment_corr,'RowNames',subject(sub_nber).flickerfreq_psd_data.label,'VariableNames',{'endogenous_freq','stim_freq'});
    end
   
end


%summarize detected endogenous frequencies:
endog_freqs=[];
num_contacts=[];
for sub_nber=1:length(subject)
    if length(subject(sub_nber).endog_peaks)<5 %means there are multiple sessions
        for s=1:length(subject(sub_nber).endog_peaks)
            for ch=1:length(subject(sub_nber).endog_peaks{s})
                if ~isempty(subject(sub_nber).endog_peaks{s}{ch})
                    endog_freqs=[endog_freqs;subject(sub_nber).endog_peaks{s}{ch}];
                    num_contacts=[num_contacts {[subject(sub_nber).subjectID ';' subject(sub_nber).flickerfreq_ssep_amp_sig{s}.Properties.RowNames{ch}]}];
                end
            end
        end
    else
        for ch=1:length(subject(sub_nber).endog_peaks)
            if ~isempty(subject(sub_nber).endog_peaks{ch})
                endog_freqs=[endog_freqs;subject(sub_nber).endog_peaks{ch}];
                num_contacts=[num_contacts {[subject(sub_nber).subjectID ';' subject(sub_nber).flickerfreq_ssep_amp_sig.Properties.RowNames{ch}]}];
            end
        end
    end
end

num_contacts=unique(num_contacts);
num_contacts=length(num_contacts);
disp(['Total number of contacts included in Figure_5D: ' num2str(num_contacts)]);

%NEED TO FIGURE OUT WHY SOME ENDOG OSCILLATIONS AMPLITUDES ARE NEGATIVE (2
%OF THEM).

figure_making('width',2,'height',2);
scatter(endog_freqs(:,1),log10(endog_freqs(:,2)),2,'k','filled','MarkerFaceAlpha',0.5);
xlabel('Frequency (Hz)');
ylabel(['Fold-change in power' newline '(log_{10} scale)']);

%save figure and associated source data:
figure_making('operation','save','filename',[figures_dir '/Figure_5D.pdf']);
writetable(cell2table(num2cell([endog_freqs(:,1),log10(endog_freqs(:,2))]),'VariableNames',{'endog_oscillation','log10_amplitude'}),[source_data '/Figure_5D.csv']);


%% identify top stimulation frequency for all channels and sessions
%including only channels of interest (i.e. ones that showed significant
%modulation to many stimulation frequencies).

correct_function=which('signal/findpeaks');
addpath(correct_function(1:end-12));

for sub_nber=1:length(subject)
    subject(sub_nber).top_modstimfreq={};
    subject(sub_nber).optimal_modstimfreq={};
    subject(sub_nber).optimal_plvstimfreq={};
    if iscell(subject(sub_nber).flickerfreq_ssep_amp_sig) %means we have multiple sessions
        for s=1:length(subject(sub_nber).flickerfreq_ssep_amp_sig)
            stim_freqs=subject(sub_nber).mod_peaks{s}.Properties.VariableNames;
            stim_freqs=str2double(regexprep(stim_freqs,'Hz-.+',''));
            [stim_freqs,stim_freq_order]=sort(stim_freqs);
            for ch=1:size(subject(sub_nber).mod_peaks{s},1)
                if sum(~isnan(subject(sub_nber).mod_peaks{s}{ch,stim_freq_order}))>nber_sigfreq_thresh %if we have x number of significant stim frequencies
                    [max_peak,temp]=max(subject(sub_nber).mod_peaks{s}{ch,stim_freq_order});
                    if isempty(temp)
                        subject(sub_nber).top_modstimfreq{s}{ch}=[];
                    else
                        subject(sub_nber).top_modstimfreq{s}{ch}=[stim_freqs(temp) max_peak];
                    end
                    
                    temp1=subject(sub_nber).mod_peaks{s}{ch,stim_freq_order};
                    temp1(isnan(temp1))=0;
                    [pks,temp]=findpeaks([0 temp1 0]); %find local maxima; add 0 at begining and end so that count maxima that are at begining and end
                    if isempty(temp)
                        subject(sub_nber).optimal_modstimfreq{s}{ch}=[];
                    else
                        subject(sub_nber).optimal_modstimfreq{s}{ch}=[stim_freqs(temp-1)' pks'];
                    end
                else
                    subject(sub_nber).top_modstimfreq{s}{ch}=[];
                    subject(sub_nber).optimal_modstimfreq{s}{ch}=[];
                end
            end
        end
    else
        stim_freqs=subject(sub_nber).mod_peaks.Properties.VariableNames;
        stim_freqs=str2double(regexprep(stim_freqs,'Hz-.+',''));
        [stim_freqs,stim_freq_order]=sort(stim_freqs);
        for ch=1:size(subject(sub_nber).mod_peaks,1)
            if sum(~isnan(subject(sub_nber).mod_peaks{ch,stim_freq_order}))>nber_sigfreq_thresh %if we have x number of significant stim frequencies
                [max_peak,temp]=max(subject(sub_nber).mod_peaks{ch,stim_freq_order});
                if isempty(temp)
                    subject(sub_nber).top_modstimfreq{ch}=[];
                else
                    subject(sub_nber).top_modstimfreq{ch}=[stim_freqs(temp) max_peak];
                end

                temp1=subject(sub_nber).mod_peaks{ch,stim_freq_order};
                temp1(isnan(temp1))=0;
                [pks,temp]=findpeaks([0 temp1 0]); %find local maxima; add 0 at begining and end so that count maxima that are at begining and end
                if isempty(temp)
                    subject(sub_nber).optimal_modstimfreq{ch}=[];
                else
                    subject(sub_nber).optimal_modstimfreq{ch}=[stim_freqs(temp-1)' pks'];
                end
            else
                subject(sub_nber).top_modstimfreq{ch}=[];
                subject(sub_nber).optimal_modstimfreq{ch}=[];
            end
        end
    end
end


%plot summary of identified optimal stim frequencies:
op_mod_stim_freqs=[];
top_mod_stim_freqs=[];
for sub_nber=1:length(subject)
    if iscell(subject(sub_nber).flickerfreq_ssep_amp_sig)
        for s=1:length(subject(sub_nber).flickerfreq_ssep_amp_sig)
            for ch=1:length(subject(sub_nber).optimal_modstimfreq{s})
                if ~isempty(subject(sub_nber).optimal_modstimfreq{s}{ch})
                    op_mod_stim_freqs=[op_mod_stim_freqs;subject(sub_nber).optimal_modstimfreq{s}{ch}];
                    top_mod_stim_freqs=[top_mod_stim_freqs;subject(sub_nber).top_modstimfreq{s}{ch}];
                end
            end
        end
    else
        for ch=1:length(subject(sub_nber).optimal_modstimfreq)
            if ~isempty(subject(sub_nber).optimal_modstimfreq{ch})
                op_mod_stim_freqs=[op_mod_stim_freqs;subject(sub_nber).optimal_modstimfreq{ch}];
                top_mod_stim_freqs=[top_mod_stim_freqs;subject(sub_nber).top_modstimfreq{ch}];
            end
        end
    end
end

figure;
scatter(top_mod_stim_freqs(:,1),log10(top_mod_stim_freqs(:,2)),'k','filled','MarkerFaceAlpha',0.5)
xlabel('Stim frequency (Hz)');
ylabel('Fold-change in power (log10 scale)');
title('Top stim frequencies');

figure;
scatter(op_mod_stim_freqs(:,1),log10(op_mod_stim_freqs(:,2)),'k','filled','MarkerFaceAlpha',0.5)
xlabel('Stim frequency (Hz)');
ylabel('Fold-change in power (log10 scale)');
title('Optimal stim frequencies');

close all;


%% identify, for each channel that has both endogenous frequency and optimal stim freq, the closest relationship between those that you can find:

plt_figure=0; %whether to plot results
plt_nber=0;
corr_matrix=[];
label_matrix={};
for sub_nber=1:length(subject)
    if iscell(subject(sub_nber).flickerfreq_ssep_amp_sig)
        for s=1:length(subject(sub_nber).flickerfreq_ssep_amp_sig)
            for ch=1:length(subject(sub_nber).endog_peaks{s})
                if ~isempty(subject(sub_nber).endog_peaks{s}{ch}) && ~isempty(subject(sub_nber).top_modstimfreq{s}{ch}) %is channel has both endog freq and top stim freq
                    [~,temp]=min(abs(subject(sub_nber).endog_peaks{s}{ch}(:,1)-subject(sub_nber).top_modstimfreq{s}{ch}(:,1)));
                    subject(sub_nber).chosen_endog_freq{s}{ch}=subject(sub_nber).endog_peaks{s}{ch}(temp,1);
                   
                    corr_matrix=[corr_matrix;subject(sub_nber).chosen_endog_freq{s}{ch} subject(sub_nber).top_modstimfreq{s}{ch}(:,1)];
                    label_matrix{end+1}=[subject(sub_nber).subjectID '; ' subject(sub_nber).flickerfreq_psd_data{s}.label{ch} ';' ...
                                            regexprep(subject(sub_nber).flickerfreq_ssep_amp_sig{s}.Properties.VariableNames{1},'.+Hz-','') ';' num2str(s)];
                    
                    if plt_figure
                        plt_nber=plt_nber+1;
                        if mod(plt_nber,6*8)==1
                            figure('Units','normalized','Position',[0 0 1 1]);
                            tiledlayout(6,8,'Padding','none','TileSpacing','none');
                        end
                        nexttile;
                        plot(subject(sub_nber).fooof_results{s}{ch}.freqs,subject(sub_nber).fooof_results{s}{ch}.power_spectrum,'k');
                        hold on;
                        plot(subject(sub_nber).fooof_results{s}{ch}.freqs,subject(sub_nber).fooof_results{s}{ch}.ap_fit,'Color',[0.5 0.5 0.5]);
                        plot(subject(sub_nber).fooof_results{s}{ch}.freqs,subject(sub_nber).fooof_results{s}{ch}.fooofed_spectrum,'r');
                        xline(subject(sub_nber).chosen_endog_freq{s}{ch},'--r');
                        
                        psd_data=subject(sub_nber).flickerfreq_psd_data{s};
                        channel=subject(sub_nber).flickerfreq_ssep_amp_sig{s}.Properties.RowNames{ch};
                        for cond=subject(sub_nber).flickerfreq_ssep_amp_sig{s}.Properties.VariableNames
                            freq_index=str2double(regexprep(cond,'Hz-.+',''));
                            freq_index=[freq_index-1 freq_index+1];
                            
                            line_color=condition_color(regexprep(cond,'.+-',''));
                            
                            frequencies=psd_data.data{strcmp(psd_data.label,channel),strcmp(psd_data.condition,cond)}{3};
                            freq_index=find(frequencies>=freq_index(1) & frequencies<=freq_index(2));
                            
                            psd_result=psd_data.data{strcmp(psd_data.label,channel),strcmp(psd_data.condition,cond)}{1};
                            plot(psd_data.data{strcmp(psd_data.label,channel),strcmp(psd_data.condition,cond)}{3}(freq_index),log10(mean(psd_result(:,freq_index))),'Color',line_color);
                            
                            x=[psd_data.data{strcmp(psd_data.label,channel),strcmp(psd_data.condition,cond)}{3}(freq_index) psd_data.data{strcmp(psd_data.label,channel),strcmp(psd_data.condition,cond)}{3}(freq_index(end:-1:1))];
                            p=patch(x,log10([mean(psd_result(:,freq_index))-std(psd_result(:,freq_index))/sqrt(size(psd_result,1)) mean(psd_result(:,freq_index(end:-1:1)))+std(psd_result(:,freq_index(end:-1:1)))/sqrt(size(psd_result,1))]),line_color,'FaceAlpha',0.2,'EdgeColor','none');
                            %p=patch(x,log10([mean(psd_result)-std(psd_result) mean(psd_result(:,end:-1:1))+std(psd_result(:,end:-1:1))]),color,'FaceAlpha',0.2,'EdgeColor','none');
                            set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                        end
                        
                        l=xline(subject(sub_nber).top_modstimfreq{s}{ch}(:,1),'--');
                        l.Color=line_color;
                        xlim([0 85.5]);
                        ylabel('Power (log_1_0)');
                    end
                                        
                end
            end
        end
    else
        for ch=1:length(subject(sub_nber).endog_peaks)
            if ~isempty(subject(sub_nber).endog_peaks{ch}) && ~isempty(subject(sub_nber).top_modstimfreq{ch}) %is channel has both endog freq and optimal stim freq
                [~,temp]=min(abs(subject(sub_nber).endog_peaks{ch}(:,1)-subject(sub_nber).top_modstimfreq{ch}(:,1)));
                subject(sub_nber).chosen_endog_freq{ch}=subject(sub_nber).endog_peaks{ch}(temp,1);
                
                corr_matrix=[corr_matrix;subject(sub_nber).chosen_endog_freq{ch} subject(sub_nber).top_modstimfreq{ch}(:,1)];
                label_matrix{end+1}=[subject(sub_nber).subjectID '; ' subject(sub_nber).flickerfreq_psd_data.label{ch} ';' ...
                    regexprep(subject(sub_nber).flickerfreq_ssep_amp_sig.Properties.VariableNames{1},'.+Hz-','')];
                                        
                if plt_figure
                    plt_nber=plt_nber+1;
                    if mod(plt_nber,6*8)==1
                        figure('Units','normalized','Position',[0 0 1 1]);
                        tiledlayout(6,8,'Padding','none','TileSpacing','none');
                    end
                    nexttile;
                    plot(subject(sub_nber).fooof_results{ch}.freqs,subject(sub_nber).fooof_results{ch}.power_spectrum,'k');
                    hold on;
                    plot(subject(sub_nber).fooof_results{ch}.freqs,subject(sub_nber).fooof_results{ch}.ap_fit,'Color',[0.5 0.5 0.5]);
                    plot(subject(sub_nber).fooof_results{ch}.freqs,subject(sub_nber).fooof_results{ch}.fooofed_spectrum,'r');
                    xline(subject(sub_nber).chosen_endog_freq{ch},'--r');

                    psd_data=subject(sub_nber).flickerfreq_psd_data;
                    channel=subject(sub_nber).flickerfreq_ssep_amp_sig.Properties.RowNames{ch};
                    for cond=subject(sub_nber).flickerfreq_ssep_amp_sig.Properties.VariableNames
                        freq_index=str2double(regexprep(cond,'Hz-.+',''));
                        freq_index=[freq_index-1 freq_index+1];

                        line_color=condition_color(regexprep(cond,'.+-',''));

                        frequencies=psd_data.data{strcmp(psd_data.label,channel),strcmp(psd_data.condition,cond)}{3};
                        freq_index=find(frequencies>=freq_index(1) & frequencies<=freq_index(2));

                        psd_result=psd_data.data{strcmp(psd_data.label,channel),strcmp(psd_data.condition,cond)}{1};
                        plot(psd_data.data{strcmp(psd_data.label,channel),strcmp(psd_data.condition,cond)}{3}(freq_index),log10(mean(psd_result(:,freq_index))),'Color',line_color);

                        x=[psd_data.data{strcmp(psd_data.label,channel),strcmp(psd_data.condition,cond)}{3}(freq_index) psd_data.data{strcmp(psd_data.label,channel),strcmp(psd_data.condition,cond)}{3}(freq_index(end:-1:1))];
                        p=patch(x,log10([mean(psd_result(:,freq_index))-std(psd_result(:,freq_index))/sqrt(size(psd_result,1)) mean(psd_result(:,freq_index(end:-1:1)))+std(psd_result(:,freq_index(end:-1:1)))/sqrt(size(psd_result,1))]),line_color,'FaceAlpha',0.2,'EdgeColor','none');
                        %p=patch(x,log10([mean(psd_result)-std(psd_result) mean(psd_result(:,end:-1:1))+std(psd_result(:,end:-1:1))]),color,'FaceAlpha',0.2,'EdgeColor','none');
                        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    end

                    l=xline(subject(sub_nber).top_modstimfreq{ch}(:,1),'--');
                    l.Color=line_color;
                    xlim([0 85.5]);
                    ylabel('Power (log_1_0)');
                end
            end
        end
    end
end

% for i=1:3
%     figure(i);
%     set(gcf,'Units','inches');
%     screenposition = get(gcf,'Position');
%     set(gcf,...
%         'PaperPosition',[0 0 screenposition(3:4)],...
%         'PaperSize',[screenposition(3:4)]);
%     print(gcf,[figures_dir '/temp_examples-entrainment' num2str(i) '.pdf'],'-dpdf','-fillpage');
%     close;
% end

figure_making('width',2,'height',2);
s=scatter(corr_matrix(:,1),corr_matrix(:,2),2,'k','filled','MarkerFaceAlpha',0.5);
hold on;
plot(linspace(2,100),linspace(2,100),'--k','Color',[0.5 0.5 0.5 0.5]);
tolerance_threshold=5;
patch([2 100 100 2],[2-tolerance_threshold 100-tolerance_threshold 100+tolerance_threshold 2+tolerance_threshold],[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',0.2);
xlim([2 100]);
ylim([2 100]);
row = dataTipTextRow('Label',label_matrix);
s.DataTipTemplate.DataTipRows(end+1) = row;
row=dataTipTextRow('dot nber',1:length(label_matrix));
s.DataTipTemplate.DataTipRows(end+1) = row;
xlabel(['Endogenous oscillation' newline 'frequency (Hz)']);
ylabel('Top stim frequency (Hz)','Color','m');

%calculate stats:
[h,p,ci,stats]=ttest(corr_matrix(:,1),corr_matrix(:,2))

disp(['p-value for ttest of the differences: ' num2str(p)]);
disp(['Percent of sig contacts with top stim freq within 5Hz of endog oscillation: ' num2str(round((sum(abs(corr_matrix(:,2)-corr_matrix(:,1))<5)/size(corr_matrix,1))*100,1)) '%']);
temp_label_matrix=arrayfun(@(x) strsplit(x{:},';'),label_matrix,'UniformOutput',false);
temp_label_matrix=arrayfun(@(x) strjoin(x{:}(1:2),';'),temp_label_matrix,'UniformOutput',false)';
disp(['Total number of contacts used in plot: ' num2str(length(unique(temp_label_matrix)))]);
temp_label_matrix=arrayfun(@(x) regexprep(x{:},';.+',''),temp_label_matrix,'UniformOutput',false);
disp(['Total number of subjects in plot: ' num2str(length(unique(temp_label_matrix)))]);
temp_label_matrix=arrayfun(@(x) strsplit(x{:},';'),label_matrix,'UniformOutput',false);
temp_label_matrix=arrayfun(@(x) strjoin(x{:}([1,3]),';'),temp_label_matrix,'UniformOutput',false)';
disp(['Total number of sessions in plot: ' num2str(length(unique(temp_label_matrix)))]);

%save figure and associated source data:
figure_making('operation','save','filename',[figures_dir '/Figure_5E.pdf']);
writetable(cell2table(num2cell(corr_matrix),'VariableNames',{'endog_oscillation_freq','top_stim_freq'}),[source_data '/Figure_5E.csv']);

%% plot illustration for all contacts:
% 
% for dot_nber=1:10%size(corr_matrix,1)
%     temp=strsplit(label_matrix{dot_nber},'; ');
%     subjectID=temp{1};
%     ch=temp{2};
%     ch_nber=find(strcmp(subject(strcmp(all_subjects,subjectID)).flickerfreq_psd_data.label,ch));
% 
%     figure;
%     subplot(2,2,1); %plot identified endogenous oscillation
%     hold on;
%     plot(subject(strcmp(all_subjects,subjectID)).fooof_results{ch_nber}.freqs,subject(strcmp(all_subjects,subjectID)).fooof_results{ch_nber}.power_spectrum);
%     plot(subject(strcmp(all_subjects,subjectID)).fooof_results{ch_nber}.freqs,subject(strcmp(all_subjects,subjectID)).fooof_results{ch_nber}.ap_fit);
%     plot(subject(strcmp(all_subjects,subjectID)).fooof_results{ch_nber}.freqs,subject(strcmp(all_subjects,subjectID)).fooof_results{ch_nber}.fooofed_spectrum);
%     xline(corr_matrix(dot_nber,1),':');
%     xlabel('Frequency (Hz)');
%     ylabel('Power (log10 scale)');
% 
%     subplot(2,2,3); %plot modulation by stim freq
%     hold on;
%     plot_ent_by_freq(subject(strcmp(all_subjects,subjectID)).flickerfreq_ssep_amp_val,subject(strcmp(all_subjects,subjectID)).flickerfreq_ssep_amp_sig,subject(strcmp(all_subjects,subjectID)).flickerfreq_ssep_amp_val.Properties.RowNames{ch_nber},1);
%     xline(corr_matrix(dot_nber,2),':');
%     xlim([0 100]);
% 
%     subplot(2,2,2); %plot PSD plot for stim freq close to endogenous frequency
%     %plot_PSD('11Hz-V',ch,subject(1).flickerfreq_psd_data,subject(1).flickerfreq_psd_data.label,subject(1).flickerfreq_psd_data.condition,[0.9100 0.4100 0.1700],1,0,1);
%     plot_PSD([num2str(corr_matrix(dot_nber,2)) 'Hz-V'],ch,subject(strcmp(all_subjects,subjectID)).flickerfreq_psd_data,subject(strcmp(all_subjects,subjectID)).flickerfreq_psd_data.label,subject(strcmp(all_subjects,subjectID)).flickerfreq_psd_data.condition,[0.9100 0.4100 0.1700],1,0,1);
%     hold on;
%     plot_PSD('Baseline',ch,subject(strcmp(all_subjects,subjectID)).flickerfreq_psd_data,subject(strcmp(all_subjects,subjectID)).flickerfreq_psd_data.label,subject(strcmp(all_subjects,subjectID)).flickerfreq_psd_data.condition,'k',1,0,1);
%     xline(corr_matrix(dot_nber,2),':');
%     xlabel('Frequency (Hz)');
%     ylabel('Power (log10 scale)');
% 
% %     subplot(2,2,4); %plot PSD plot for stim freq close to endogenous frequency
% %     plot_PSD('57Hz-V',ch,subject(1).flickerfreq_psd_data,subject(1).flickerfreq_psd_data.label,subject(1).flickerfreq_psd_data.condition,[0.9100 0.4100 0.1700],1,0,1);
% %     hold on;
% %     plot_PSD('Baseline',ch,subject(1).flickerfreq_psd_data,subject(1).flickerfreq_psd_data.label,subject(1).flickerfreq_psd_data.condition,'k',1,0,1);
% %     xline(57,':');
% %     xlabel('Frequency (Hz)');
% %     ylabel('Power (log10 scale)');
% end

%% plot examples:

examples=readtable([figures_dir '/NatureComm2024_examples.xlsx']);
examples=table2array(examples(startsWith(examples.example,'Figure_5'),:));
examples=examples(:,[2 5 4 1]);
[~,temp_index]=sort(examples(:,4));
examples=examples(temp_index,:);

figure_making('width',4.5,'height',2);
tiledlayout(2,3,'TileSpacing','normal','Padding','normal');

for i=1:size(examples,1)
    
    nexttile(i);
    
    sub_nber=strcmp(all_subjects,examples{i,1});
    if isempty(examples{i,3})
        ch_nber=strcmp(subject(sub_nber).flickerfreq_psd_data.label,examples{i,2});
        modality=regexprep(subject(sub_nber).flickerfreq_ssep_amp_sig.Properties.VariableNames{1},'.+-','');
        conditions=subject(sub_nber).flickerfreq_ssep_amp_sig.Properties.VariableNames;
        PSD_data_obj=subject(sub_nber).flickerfreq_psd_data;
        plot(subject(sub_nber).fooof_results{ch_nber}.freqs,subject(sub_nber).fooof_results{ch_nber}.ap_fit,'LineStyle','-','Color',[0.5 0.5 0.5]);
        hold on;
        %plot(subject(sub_nber).fooof_results{ch_nber}.freqs,subject(sub_nber).fooof_results{ch_nber}.fooofed_spectrum,'LineStyle','-','Color','m');
    elseif ~isempty(examples{i,3})
        ch_nber=strcmp(subject(sub_nber).flickerfreq_psd_data{str2num(examples{i,3})}.label,examples{i,2});
        modality=regexprep(subject(sub_nber).flickerfreq_ssep_amp_sig{str2num(examples{i,3})}.Properties.VariableNames{1},'.+-','');
        conditions=subject(sub_nber).flickerfreq_ssep_amp_sig{str2num(examples{i,3})}.Properties.VariableNames;
        PSD_data_obj=subject(sub_nber).flickerfreq_psd_data{str2num(examples{i,3})};
        plot(subject(sub_nber).fooof_results{str2num(examples{i,3})}{ch_nber}.freqs,subject(sub_nber).fooof_results{str2num(examples{i,3})}{ch_nber}.ap_fit,'LineStyle','-','Color',[0.5 0.5 0.5]);
        hold on;
        %plot(subject(sub_nber).fooof_results{str2num(examples{i,3})}{ch_nber}.freqs,subject(sub_nber).fooof_results{str2num(examples{i,3})}{ch_nber}.fooofed_spectrum,'LineStyle','-','Color','m');
    end
    
    xlim([0 85.5]);
    line_color=condition_color(modality);
    for cond=conditions
        freq_index=str2double(regexprep(cond,'Hz-.+',''));
        freq_index=[freq_index-1 freq_index+1];

        frequencies=PSD_data_obj.data{strcmp(PSD_data_obj.label,examples{i,2}),strcmp(PSD_data_obj.condition,cond)}{3};
        freq_index=find(frequencies>=freq_index(1) & frequencies<=freq_index(2));

        psd_result=PSD_data_obj.data{strcmp(PSD_data_obj.label,examples{i,2}),strcmp(PSD_data_obj.condition,cond)}{1};
        plot(PSD_data_obj.data{strcmp(PSD_data_obj.label,examples{i,2}),strcmp(PSD_data_obj.condition,cond)}{3}(freq_index),log10(mean(psd_result(:,freq_index))),'Color',line_color);

        x=[PSD_data_obj.data{strcmp(PSD_data_obj.label,examples{i,2}),strcmp(PSD_data_obj.condition,cond)}{3}(freq_index) PSD_data_obj.data{strcmp(PSD_data_obj.label,examples{i,2}),strcmp(PSD_data_obj.condition,cond)}{3}(freq_index(end:-1:1))];
        p=patch(x,log10([mean(psd_result(:,freq_index))-std(psd_result(:,freq_index))/sqrt(size(psd_result,1)) mean(psd_result(:,freq_index(end:-1:1)))+std(psd_result(:,freq_index(end:-1:1)))/sqrt(size(psd_result,1))]),line_color,'FaceAlpha',0.2,'EdgeColor','none');
        %p=patch(x,log10([mean(psd_result)-std(psd_result) mean(psd_result(:,end:-1:1))+std(psd_result(:,end:-1:1))]),color,'FaceAlpha',0.2,'EdgeColor','none');
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

    end
    
    plot_PSD('Baseline',examples{i,2},PSD_data_obj,PSD_data_obj.label,PSD_data_obj.condition,'k',1,0,1);
    
    if isempty(examples{i,3})
        xline(subject(sub_nber).top_modstimfreq{ch_nber}(1,1),'--','Color',[line_color 0.5]);
        xline(subject(sub_nber).chosen_endog_freq{ch_nber},'--','Color','k');
    elseif ~isempty(examples{i,3})
        xline(subject(sub_nber).top_modstimfreq{str2num(examples{i,3})}{ch_nber}(1,1),'--','Color',[line_color 0.5]);
        xline(subject(sub_nber).chosen_endog_freq{str2num(examples{i,3})}{ch_nber},'--','Color','k');
    end
    
    if i==1
        ylabel(['Power' newline '(log_{10} scale)']);
    elseif i==3
        %add legend:
        ax=gca;
        text(70,ax.YLim(2)-range(ax.YLim)/2.5,['\color[rgb]{' num2str(condition_color('V')) '}V' newline '\color[rgb]{' num2str(condition_color('A')) '}A' newline '\color[rgb]{0 0 0}BL' newline '\color[rgb]{0.5 0.5 0.5}1/f']);
    end
    box off;
    
    %plot ent_by_freq and plv response:
    nexttile(i+3);
    yyaxis left;
    if isempty(examples{i,3})
        plot_ent_by_freq(subject(strcmp(all_subjects,examples{i,1})).flickerfreq_ssep_amp_val,subject(strcmp(all_subjects,examples{i,1})).flickerfreq_ssep_amp_sig,examples{i,2},1);
    elseif ~isempty(examples{i,3})
        plot_ent_by_freq(subject(sub_nber).flickerfreq_ssep_amp_val{str2num(examples{i,3})},subject(sub_nber).flickerfreq_ssep_amp_sig{str2num(examples{i,3})},examples{i,2},1);
    end
    xlabel([]);
    if i==1
        ylabel(['Fold change' newline 'in power']);
    else
        ylabel([]);
    end
    
    yyaxis right;
    if isempty(examples{i,3})
        plot_ent_by_freq(subject(sub_nber).flickerfreq_ssep_plv_amp,subject(strcmp(all_subjects,examples{i,1})).flickerfreq_ssep_plv_sig,examples{i,2},1,'plv');
    elseif ~isempty(examples{i,3})
        plot_ent_by_freq(subject(sub_nber).flickerfreq_ssep_plv_amp{str2num(examples{i,3})},subject(sub_nber).flickerfreq_ssep_plv_sig{str2num(examples{i,3})},examples{i,2},1,'plv');
    end
    xlabel([]);
    if i==3
        ylabel('PLV');
    else
        ylabel([]);
    end
    
    if isempty(examples{i,3})
        xline(subject(sub_nber).top_modstimfreq{ch_nber}(1,1),'--','Color',[line_color 0.5]);
        %xline(subject(sub_nber).chosen_endog_freq{ch_nber},'--','Color','m');
    elseif ~isempty(examples{i,3})
        xline(subject(sub_nber).top_modstimfreq{str2num(examples{i,3})}{ch_nber}(1,1),'--','Color',[line_color 0.5]);
        %xline(subject(sub_nber).chosen_endog_freq{str2num(examples{i,3})}{ch_nber},'--','Color','m');
    end
    
    title([]);
    box off;
    ax=gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    line([ax.XLim(2) ax.XLim(2)],ax.YLim,'color','w','LineWidth',1);
    line([ax.XLim(2) ax.XLim(2)],ax.YLim,'color','k','LineWidth',0.5,'LineStyle','--');
end
han=axes(gcf,'visible','off'); 
han.XLabel.Visible='on';
xlabel(han,'Stimulation frequency (Hz)');

figure_making('operation','save','filename',[figures_dir '/Figure_5B.pdf']);

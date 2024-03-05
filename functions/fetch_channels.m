%Gets list of clinical channel labels, and corresponding montage channel
%numbers in some cases. Additionally, has mechanism to proofcheck that have
%correct channels labels in various files.
%
%fnames: structure that may contain:
%*root_dir
%*subjectID
%*task
%*ses
%2024/02/25

function varargout=fetch_channels(fnames,request)
    %fetch implant list channel labels
    all_channel_labels=readtable([fnames.root_dir '/stg-raw/sub-' fnames.subjectID '/anat/sub-' fnames.subjectID '_stg-raw_desc-implantList_nat-list.csv']);
    %check for whitespaces in implant list labels
    if any(arrayfun(@(x) any(isspace(x{:})),all_channel_labels.label))
        error(['There are whitespaces in the following labels of implantList: ' strjoin(all_channel_labels.label(arrayfun(@(x) any(isspace(x{:})),all_channel_labels.label)),', ')]);
    end
    
    %perform various tasks:
    if strcmp(request,'all') %fetch all electrode contacts
        output_arg=1;
        channels=all_channel_labels;
    elseif strcmp(request,'all_rec') %fetch all clinical contacts that are in montage
        output_arg=1;
        channels=readtable([fnames.root_dir '/stg-raw/sub-' fnames.subjectID '/anat/sub-' fnames.subjectID '_stg-raw_desc-montage_nat-list.csv']);
        if ~all(ismember(channels.label,table2array(all_channel_labels))) %check that all montage labels are included in list of potential channel labels
            error(['The following montage labels do not exist in list of all clinical channel labels: ' strjoin(channels.label(~ismember(channels.label,table2array(all_channel_labels))),', ')]);
        end
    elseif strcmp(request,'-chExclude') %do not include excluded channels
        output_arg=1;
        channels=readtable([fnames.root_dir '/stg-raw/sub-' fnames.subjectID '/anat/sub-' fnames.subjectID '_stg-raw_desc-montage_nat-list.csv']);
        if ~all(ismember(channels.label,table2array(all_channel_labels))) %check that all montage labels are included in list of potential channel labels
            error(['The following montage labels do not exist in list of all clinical channel labels: ' strjoin(channels.label(~ismember(channels.label,table2array(all_channel_labels))),', ')]);
        end
        chExclude_labels=readcell([fnames.root_dir '/stg-raw/sub-' fnames.subjectID '/anat/sub-' fnames.subjectID '_stg-raw_desc-chExclude_nat-list.txt'],'CommentStyle','#');
        if ~all(ismember(chExclude_labels,table2array(all_channel_labels))) %check that all chExclude labels are included in list of potential channel labels
            error(['The following chExclude labels do not exist in list of all clinical channel labels: ' strjoin(chExclude_labels(~ismember(chExclude_labels,table2array(all_channel_labels))),', ')]);
        else
            channels=channels(~ismember(channels.label,chExclude_labels),:);
        end
    elseif strcmp(request,'-chNoisy') %do no include noisy channels
        output_arg=1;
        channels=fetch_channels(fnames,'-chExclude'); %get list of channels without the exclude channels
        chNoisy_labels=readcell([fnames.root_dir '/stg-preproc/sub-' fnames.subjectID '/task-' fnames.task '/ses-' fnames.ses '/sub-' fnames.subjectID '_stg-preproc_task-' fnames.task '_ses-' fnames.ses '_desc-noisy-ch_nat-list.csv'],'FileType','text');
        if ~all(ismember(chNoisy_labels,table2array(all_channel_labels))) %check that all chNoisy labels are included in list of potential channel labels
            error(['The following chNoisy labels do not exist in list of all clinical channel labels: ' strjoin(chNoisy_labels(~ismember(chNoisy_labels,table2array(all_channel_labels))),', ')]);
        else
            channels=channels(~ismember(channels.label,chNoisy_labels),:);
        end
    elseif strcmp(request,'proofcheck-voxtool') %check that labels outputted by voxtool are correct
        output_arg=0;
        voxtool_labels=extract_vt_json([fnames.root_dir '/stg-preproc/sub-' fnames.subjectID '/anat/voxtool-output/sub-' fnames.subjectID '_e-vox-coords_ct.json']);
        voxtool_labels=voxtool_labels.label;
        if ~all(ismember(voxtool_labels,table2array(all_channel_labels)))
            error(['The following voxtool labels are wrong: ' strjoin(voxtool_labels(~ismember(voxtool_labels,table2array(all_channel_labels))),', ') '.']);
        end
    elseif strcmp(request,'proofcheck-clinicallabels') %check that neurologist anatomical assignment list  has correct labels
        output_arg=0;
        clinical_labels=readtable([fnames.root_dir '/stg-raw/sub-' fnames.subjectID '/anat/sub-' fnames.subjectID '_clinicalElectrodeLocalization.csv']);
        clinical_labels=clinical_labels.label;
        if ~all(ismember(clinical_labels,table2array(all_channel_labels)))
            error(['The following clinical anatomical assignment list labels are wrong: ' strjoin(clinical_labels(~ismember(clinical_labels,table2array(all_channel_labels))),', ') '.']);
        end
    elseif strcmp(request,'proofcheck-EDFlabels') %compare montage labels to EDF file labels
        output_arg=0;
        channels=cell([256,3]); %initialize comparison table
        channels(:,1)=num2cell((1:256)'); %place channel numbers
        montage_channels=readtable([fnames.root_dir '/stg-raw/sub-' fnames.subjectID '/anat/sub-' fnames.subjectID '_stg-raw_desc-montage_nat-list.csv']); %fetch clinical montage
        if ~all(ismember(montage_channels.label,table2array(all_channel_labels))) %check that all montage labels are included in list of potential channel labels
            error(['The following montage labels do not exist in list of all clinical channel labels: ' strjoin(montage_channels.label(~ismember(montage_channels.label,table2array(all_channel_labels))),', ')]);
        end
        channels(montage_channels.channel,2)=montage_channels.label; %add montage labels
        edf_labels=ft_read_header([fnames.root_dir '/stg-raw/sub-' fnames.subjectID '/task-' fnames.task '/ses-' fnames.ses '/sub-' fnames.subjectID '_stg-raw_task-' fnames.task '_ses-' fnames.ses '_nat-ieeg.EDF']); %fetch edf file labels
        if size(edf_labels.label,1)<256 %means only 128 channels were used in clinical amplifier (based on EDF file)
            if all(arrayfun(@(x) isempty(x{:}),channels(129:end,2))) %double-check that none of the channels from 129-256 have labels in the montage
                channels(129:end,:)=[]; %remove those extra channels
                edf_labels=edf_labels.label(1:128);
            else
                error('Some channels in montage have labels, but were apparently not recorded from in EDF file...');
            end
        else %means 256 channels were used in clinical amplifier
            edf_labels=edf_labels.label(1:256);
        end
        channels(:,3)=edf_labels; %add edf labels
        channels=cell2table(channels,'VariableNames',{'channel','montage_label','edf_label'});
        channels(arrayfun(@(x) isempty(x{:}),channels.montage_label),:)=[]; %remove channels that are not part of the montage
        temp=table2array(channels(:,{'montage_label','edf_label'}));
        channels(strcmp(temp(:,1),temp(:,2)),:)=[]; %remove channels for which montage and edf labels match
        if ~isempty(channels) %if we have channels for which montage and edf labels don't match, warn user and return those channels
            warning('Some edf labels do not match montage labels:');
            disp(channels);
        end
    elseif startsWith(request,'patho_channels') %fetch list of channels associated with some pathology (i.e. abnormal tissue, SOZ or interictal spikes)
        output_arg=1;
        channels=cell([0 2]);
        fileID=fopen([fnames.root_dir '/stg-raw/sub-' fnames.subjectID '/anat/sub-' fnames.subjectID '_stg-raw_desc-chEpileptiform_nat-list.txt'],'r');
        number_hashtags=0; %assumes that channels are in that order: abnormal tissue, SOZ, IEDs
        while ~feof(fileID) %while we have not reached end of file
            current_line=fgetl(fileID); %get next line
            if startsWith(current_line,'#')
                number_hashtags=number_hashtags+1;
                switch number_hashtags %channels are ordered by those implicated in abnormal tissue first, then soz, then ieds
                    case 1
                        feature_name='abnormal';
                    case 2
                        feature_name='soz';
                    case 3
                        feature_name='ied';
                end
            elseif isempty(current_line)
                %do nothing
            else %means we have a label
                channels(end+1,:)={current_line,feature_name};
            end
        end
        fclose(fileID);
        
        channels=array2table(channels,'VariableNames',{'label','feature'});
        
        if ~all(ismember(channels.label,table2array(all_channel_labels))) %check that all retrieved labels are included in list of potential channel labels
            error(['The following retrieved labels do not exist in list of all clinical channel labels: ' strjoin(channels.label(~ismember(channels.label,table2array(all_channel_labels))),', ')]);
        end
    end
    
    if output_arg==1
        varargout{1}=channels;
    end
end

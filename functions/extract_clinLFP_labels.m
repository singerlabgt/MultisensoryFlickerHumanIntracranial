%Organizes electrode contact labels by depth electrode.
%05/27/2021: check that this works for new Emory electrode nomenclature
%(i.e. "A'1" for example).
%2024/02/25

function clinLFP_contacts=extract_clinLFP_labels(labels_array)
    %remove micro contacts if there are any:
    labels_array(contains(labels_array,'Micro'))=[];
    
    %reorder electrode labels if needed:
    [labels_array,~]=reorder_elabels(labels_array);
    
    %get electrode names:
    temp_labels=cellfun(@(x) strsplit(x,'-'),labels_array,'UniformOutput',false); %first do this in case we have a bipolar or laplacean montage
    if length(temp_labels{1})>1 && ~any(cell2mat(cellfun(@(x) isempty(regexp(x{1,1}(end),'\d')),temp_labels,'UniformOutput',false))) %added 2nd condition because on very rare cases, can get weird contact labels with '-', such as 1Ld-2Rd1 (in non-Laplacian or bipolar montage)
        temp_labels=cellfun(@(x) x{1,:},temp_labels,'UniformOutput',false)';
    else
        temp_labels=labels_array;
    end
    
    temp=regexp(temp_labels,'\d*'); %get index of characters where a number starts
    label_root={};
    for i=1:length(temp_labels)
        label_root(end+1,:)={temp_labels{i}(1:temp{i}(end)-1) temp_labels{i}(temp{i}(end):end)};
    end
    
    electrodes_labels=unique(label_root(:,1),'stable');
    
    %place all electrode labels in an organized structure:
    clinLFP_contacts=struct();
    for i=1:length(electrodes_labels)
        clinLFP_contacts(i).depth_electrode_name=electrodes_labels{i};
        current_contacts=labels_array(ismember(label_root,clinLFP_contacts(i).depth_electrode_name));
        clinLFP_contacts(i).channel_names=current_contacts;
    end
end

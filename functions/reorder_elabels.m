%Reorders electrode contact labels.
%05/27/2021: checked that this works for new Emory electrode nomenclature (i.e. "A'1"
%for example).
%2024/02/26

function [labels_array,labels_index]=reorder_elabels(labels_array)
    
    temp_labels=cellfun(@(x) strsplit(x,'-'),labels_array,'UniformOutput',false); %first do this in case we have a bipolar or laplacean montage
    if length(temp_labels{1})>1 && ~isempty(regexp(temp_labels{1}{1}(end),'\d')) %added 2nd condition because on very rare cases, can get weird contact labels with '-', such as 1Ld-2Rd1 (in non-Laplacian or bipolar montage)
        temp_labels=cellfun(@(x) x{1,:},temp_labels,'UniformOutput',false)';
    else
        temp_labels=labels_array;
    end

    %reorder electrode labels:
    temp=regexp(temp_labels,'\d*'); %get index of characters where a number starts
    numbers={};
    for i=1:length(temp_labels)
        numbers(end+1,:)={temp_labels{i}(1:temp{i}(end)-1) temp_labels{i}(temp{i}(end):end)};
    end
    
    %need this to order labels:
    numbers=[repmat({'0'},[size(numbers,1),1]) numbers(:,1) repmat({'0'},[size(numbers,1),1]) numbers(:,2)];
    
    %in case we have micro contacts:
    if any(contains(numbers(:,2),'Micro')) %means we have micro wires
        numbers(contains(numbers(:,2),'Micro'),3)={'1'};
        numbers(:,2)=regexprep(numbers(:,2),'Micro','');
    end
    
    %if there are electrode that have numbers in begining of their labels,
    %reorder them in increasing number:
    for i=1:length(numbers(:,2))
        temp=regexp(numbers(i,2),'\d*','match');
        if length(temp{:})>1 %added this because on very rare cases, can get weird contact labels with '-', such as 1Ld-2Rd1 (in non-Laplacian or bipolar montage)
            temp={temp{:}(1)};
        end
        if ~isempty(temp{:}) && startsWith(numbers(i,2),temp{1,1})
            numbers(i,1)=temp{1,1};
        end
    end
    
    %reorder labels:
    numbers(:,[1 3:4])=cellfun(@(x) str2num(x),numbers(:,[1 3:4]),'UniformOutput',false);
    [~,labels_index]=sortrows(numbers,1:size(numbers,2));
    labels_array=labels_array(labels_index);
    
end
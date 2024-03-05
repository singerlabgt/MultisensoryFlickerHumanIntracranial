%Fetches unique subject IDs for various experiments.

function [subjectIDs,sessions]=fetch_flicker_subjectIDs(root_dir,task_name)
    metadata_tbl=readtable([root_dir '/FlickerStudyMetadata.xlsx'],'Sheet','Sessions','PreserveVariableNames',1);
    metadata_tbl=metadata_tbl(:,{'Subject_ID','Experiment','Session','Session_version','Include in analysis?'});
    
    if strcmp(task_name,'all') %gives subject IDs for subjects who completed at least 1 task
        sessions=metadata_tbl(logical(metadata_tbl.("Include in analysis?")),{'Subject_ID','Experiment','Session','Session_version'});
    else
        sessions=metadata_tbl(all([logical(metadata_tbl.("Include in analysis?")) strcmp(metadata_tbl.Experiment,task_name)]'),{'Subject_ID','Experiment','Session','Session_version'});
    end
    
    %rename variable names:
    sessions.Properties.VariableNames={'sub','task','ses','modality'};
    
    %convert session numbers to proper strings:
    sessions.ses=cellstr(strcat('0',num2str(sessions.ses)));
    
    %get unique subject IDs for this:
    subjectIDs=unique(sessions.sub);
    
end

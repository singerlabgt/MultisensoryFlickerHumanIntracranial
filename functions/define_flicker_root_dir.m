%Defines where flicker data is, as well as where dependent repositories
%are.
%2024/02/25

function [flicker_root_dir, repo]=define_flicker_root_dir()
    
    %root directories:
    if ispc
        root_dir='Z:';
        repo_dir='C:/Users/lblanpa/Documents/Code';
    elseif ismac
        root_dir='/Volumes/groups/';
        repo_dir='/Users/lou.blanpain/Documents/Coding/code/';
    end

    %define where flicker data is:
    flicker_root_dir=[root_dir '/WillieFlickerStudy/NatureComm2024'];
    
    %define where dependent repositories are:
    repo_dir=[repo_dir '/WillieLabRepositories'];
    repo=struct;
    repo.fieldtrip=[repo_dir '/Analysis_Fieldtrip'];
    repo.chronux=[repo_dir '/Analysis_Chronux'];
    repo.fooof=[repo_dir '/../fooof_mat-main'];
    repo.circstats=[repo_dir '/Analysis_CircStats'];
    repo.exportfig=[repo_dir '/export_fig'];
    repo.violinplot=[repo_dir '/Violinplot-Matlab'];
    repo.venn=[repo_dir '/Analysis_CommonCode/downloadedFunctions/venn'];
    
end

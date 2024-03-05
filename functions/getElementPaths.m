
% Gets a list of strings specifying paths of elements in inputFolder that
% match pattern.

function filePathsOfInterest=getElementPaths(inputFolder,pattern,varargin)
    
    if isempty(varargin)
        typeOfElement='file';
    elseif length(varargin)==1
        typeOfElement=varargin{1};
    end
    
    elements_of_interest=dir(inputFolder);
    
    if strcmp(typeOfElement,'folder')
        elements_of_interest=elements_of_interest([elements_of_interest.isdir]); %only keep directories
        elements_of_interest=elements_of_interest(~contains({elements_of_interest.name},'.')); %remove . and ..
    elseif strcmp(typeOfElement,'file')
        elements_of_interest=elements_of_interest(~[elements_of_interest.isdir]); %only keep directories
    end
    
    if contains(pattern,'*')
        pattern=regexprep(pattern,'*','.+');
    end
    
    strings_containing_patterns=regexp({elements_of_interest.name},pattern); %find which strings contain the pattern
    strings_containing_patterns=cellfun(@(x) length(x),strings_containing_patterns); %find how many matches were found, and throw error if more than one matches were found
    if any(strings_containing_patterns>1)
        error('Set of patterns is matched more than once in some of the strings');
    end
    strings_containing_patterns=logical(strings_containing_patterns);
    
    elements_of_interest=elements_of_interest(strings_containing_patterns); %keep only keep elements that fit name pattern
    
    filePathsOfInterest=string(strcat(inputFolder,filesep,{elements_of_interest.name}))';
    
end
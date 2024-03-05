%Fetches mni coordinates for all contacts across subjects that are within
%subregions of interest (contrast of 2 broader regions), with their flicker response amplitude and
%significance.
%Requires structure subject to have fields .anat, .flickerneuro_ssep_amp_val
%and .flickerneuro_ssep_amp_sig.

function [mat1,mat2,mat_labels1,mat_labels2]=fetch_anat_response(subject,anat_labels1,anat_labels2,condition)

    mat1=[];
    mat2=[];
    mat_labels1={};
    mat_labels2={};

    for i=1:length(subject) %for each subject
        for j=1:length(subject(i).anat.electrodes_info.labels) %for each electrode
            if any(startsWith(subject(i).flickerneuro_ssep_amp_val.Properties.RowNames,[subject(i).anat.electrodes_info.labels{j} '-'])) %if we actually ran that channel through our analysis
                if endsWith(subject(i).anat.electrodes_info.anatlabels{j,'fs_aparcaseg'},anat_labels1) %means contact is in visual region (according to freesurfer)
                    mat1(end+1,:)=[subject(i).anat.fs_coords_mni(j,:) subject(i).flickerneuro_ssep_amp_val{startsWith(subject(i).flickerneuro_ssep_amp_val.Properties.RowNames,[subject(i).anat.electrodes_info.labels{j} '-']),condition} subject(i).flickerneuro_ssep_amp_sig{startsWith(subject(i).flickerneuro_ssep_amp_sig.Properties.RowNames,[subject(i).anat.electrodes_info.labels{j} '-']),condition}]; %add mni position and flicker amplitude to visual_regions_mat
                    mat_labels1{end+1}=[subject(i).subjectID ';' subject(i).anat.electrodes_info.labels{j}];
                elseif endsWith(subject(i).anat.electrodes_info.anatlabels{j,'fs_aparcaseg'},anat_labels2) %means contact is in frontal lobe (according to freesurfer)
                    mat2(end+1,:)=[subject(i).anat.fs_coords_mni(j,:) subject(i).flickerneuro_ssep_amp_val{startsWith(subject(i).flickerneuro_ssep_amp_val.Properties.RowNames,[subject(i).anat.electrodes_info.labels{j} '-']),condition} subject(i).flickerneuro_ssep_amp_sig{startsWith(subject(i).flickerneuro_ssep_amp_sig.Properties.RowNames,[subject(i).anat.electrodes_info.labels{j} '-']),condition}]; %add mni position and zscore to PFC mat
                    mat_labels2{end+1}=[subject(i).subjectID ';' subject(i).anat.electrodes_info.labels{j}];
                end
            end
        end
    end
    
end
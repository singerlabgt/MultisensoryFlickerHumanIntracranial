%Calculates significance value using repeated shuffling of trials.
%2024/02/25

function p_value=pval_randomshuffle(two_col_vector,num_shuffles)
    x=0;
    for i=1:num_shuffles %for number of shufflings
        %randomize trials conditions
        temp=two_col_vector(:);
        temp=temp(randperm(length(temp)));
        temp=[temp(1:size(two_col_vector,1)) temp(size(two_col_vector,1)+1:end)];
        
        %assess whether this shuffle iteration has a mean value greater than measured:
        if mean(temp(:,1))-mean(temp(:,2))>mean(two_col_vector(:,1))-mean(two_col_vector(:,2))
            x=x+1;
        end
    end
    
    %estimate p-value:
    p_value=x/num_shuffles;
end

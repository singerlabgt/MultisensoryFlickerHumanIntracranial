%Segment flicker trials into 2-pulse trials.
%2024/02/26

function [conditionIndex,result]=converttimestamps3(trials,field,numberCycles,repeatCycles)
    %create index of conditions:
    conditionIndex={'5.5Hz-Baseline' '40Hz-Baseline' '80Hz-Baseline' ...
            '5.5Hz-RV' '5.5Hz-RAV' '5.5Hz-RA' ...
            '40Hz-RV' '40Hz-RAV' '40Hz-RA' ...
            '80Hz-RV' '80Hz-RAV' '80Hz-RA' ...
            '5.5Hz-V' '5.5Hz-AV' '5.5Hz-A' ...
            '40Hz-V' '40Hz-AV' '40Hz-A' ...
            '80Hz-V' '80Hz-AV' '80Hz-A'};
        
    %add 40Hz-AV_occluded condition index if that condition exists:
    if any(strcmp(trials.condition_code,'40Hz-AV_occluded')) 
        conditionIndex{end+1}='40Hz-AV_occluded';
    end
    
    %add 40Hz-V_occluded and 40Hz-A_occluded condition indices if those conditions exist
    if any(strcmp(trials.condition_code,'occluded_40Hz-V'))
        conditionIndex{end+1}='occluded_40Hz-V';
        conditionIndex{end+1}='occluded_40Hz-A';
    end
    
    %segment into 2-pulse trials:
    result={};
    pulse_being_processed=1;
    currentConditionIndex=[];
    conditionIndexRow=num2cell(ones(1,length(conditionIndex)));
    while pulse_being_processed<=size(trials.reconstructedPulses_identities,1) %for each pulse we need to process
        currentConditionNumber=trials.reconstructedPulses_identities(pulse_being_processed,1); %find condition code for that pulse
        currentCondition=trials.condition_code(currentConditionNumber); %find corresponding condition name
        if ~contains(currentCondition,"R") && ~contains(currentCondition,"Baseline") && ~contains(currentCondition,"occluded") %if this is not a condition involving random, baseline or occluded
            currentConditionIndex=find(strcmp(conditionIndex,currentCondition)); %find current condition index
            %find indices of all pulses for this trial:
            for i=pulse_being_processed:size(trials.reconstructedPulses_identities,1)
                if trials.reconstructedPulses_identities(i,1)~=trials.reconstructedPulses_identities(pulse_being_processed)
                    break; %break once we reach a pulse belonging to a new condition (means we are moving on to next trial)
                end
            end
            block_indices=pulse_being_processed:i-1; %gives indices of pulses belonging to current trial
            
            %segment trial:
            if ~repeatCycles %if not repeating data from cycles
                if mod(length(block_indices),numberCycles)==0 %if number of pulses is a multiple of number of cycles we want to represent
                    for i=block_indices(1):numberCycles:block_indices(end) %NEED TO FIX SO THAT DON'T IGNORE THE LAST FEW CYCLES IN A TRIAL 
                        if i~=block_indices(end)-numberCycles+1
                            result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[trials.(field).reconstructedPulses_timestamps(i,1) trials.(field).reconstructedPulses_timestamps(i+numberCycles,1)-1 trials.reconstructedPulses_identities(pulse_being_processed,:)];
                            conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                        elseif i==block_indices(end)-numberCycles+1 %means we're on the last iteration of i- we need to estimate the end of the off period of the last pulse
                            temp=str2double(regexprep(currentCondition,'Hz.+',''));
                            result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[trials.(field).reconstructedPulses_timestamps(i,1) trials.(field).reconstructedPulses_timestamps(i+numberCycles-1,2)+round(trials.(field).sampleRate/temp/2) trials.reconstructedPulses_identities(pulse_being_processed,:)];
                            conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                        end
                    end
                elseif mod(length(block_indices),numberCycles)~=0 %if number of pulses is not a multiple of number of cycles we want to represent
                    for i=block_indices(1):numberCycles:block_indices(end)-mod(length(block_indices),numberCycles) %in this case, we'll ignore the last pulses not included in the the multiple of number of cycles we want
                        result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[trials.(field).reconstructedPulses_timestamps(i,1) trials.(field).reconstructedPulses_timestamps(i+numberCycles,1)-1 trials.reconstructedPulses_identities(pulse_being_processed,:)];
                        conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                    end
                end
            elseif repeatCycles %if repeating data from cycles
                for i=block_indices(1:end-(numberCycles-1)) 
                    if i~=block_indices(end-(numberCycles-1)) 
                        result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[trials.(field).reconstructedPulses_timestamps(i,1) trials.(field).reconstructedPulses_timestamps(i+numberCycles,1)-1 trials.reconstructedPulses_identities(pulse_being_processed,:)];
                        conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                    elseif i==block_indices(end-(numberCycles-1))  %if we're on last iteration of i, we need to estimate the end of the off period of the last pulse
                        temp=str2double(regexprep(currentCondition,'Hz.+',''));
                        result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[trials.(field).reconstructedPulses_timestamps(i,1) trials.(field).reconstructedPulses_timestamps(i+numberCycles-1,2)+round(trials.(field).sampleRate/temp/2) trials.reconstructedPulses_identities(pulse_being_processed,:)];
                        conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                    end
                end
            end
            pulse_being_processed=pulse_being_processed+length(block_indices);
        elseif strcmp("Baseline",currentCondition) %if current condition is baseline %NEED TO PRROF-CHECK/ WORK ON BASELINE AND RANDOM
            for condition={'5.5Hz-Baseline' '40Hz-Baseline' '80Hz-Baseline'}
                 switch char(condition)
                    case '5.5Hz-Baseline'
                        if ~repeatCycles
                            trialLength=(1/5.5)*trials.(field).sampleRate*numberCycles;
                            temp=trials.(field).reconstructedPulses_timestamps(pulse_being_processed,1):trialLength:trials.(field).reconstructedPulses_timestamps(pulse_being_processed,2);
                            currentConditionIndex=find(strcmp(conditionIndex,condition));
                            for j=1:length(temp)
                                result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[temp(j) temp(j)+trialLength-1 trials.reconstructedPulses_identities(pulse_being_processed,:)];
                                conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                            end
                        elseif repeatCycles
                            trialLength=(1/5.5)*trials.(field).sampleRate*numberCycles;
                            temp=trials.(field).reconstructedPulses_timestamps(pulse_being_processed,1):trialLength/numberCycles:trials.(field).reconstructedPulses_timestamps(pulse_being_processed,2);
                            currentConditionIndex=find(strcmp(conditionIndex,condition));
                            for j=1:length(temp)-1
                                result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[temp(j) temp(j)+trialLength-1 trials.reconstructedPulses_identities(pulse_being_processed,:)];
                                conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                            end
                        end
                    case '40Hz-Baseline'
                        if ~repeatCycles
                            trialLength=(1/40)*trials.(field).sampleRate*numberCycles;
                            temp=trials.(field).reconstructedPulses_timestamps(pulse_being_processed,1):trialLength:trials.(field).reconstructedPulses_timestamps(pulse_being_processed,2);
                            currentConditionIndex=find(strcmp(conditionIndex,condition));
                            for j=1:length(temp)
                                result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[temp(j) temp(j)+trialLength-1 trials.reconstructedPulses_identities(pulse_being_processed,:)];
                                conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                            end
                        elseif repeatCycles
                            trialLength=(1/40)*trials.(field).sampleRate*numberCycles;
                            temp=trials.(field).reconstructedPulses_timestamps(pulse_being_processed,1):trialLength/numberCycles:trials.(field).reconstructedPulses_timestamps(pulse_being_processed,2);
                            currentConditionIndex=find(strcmp(conditionIndex,condition));
                            for j=1:length(temp)-1
                                result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[temp(j) temp(j)+trialLength-1 trials.reconstructedPulses_identities(pulse_being_processed,:)];
                                conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                            end
                        end
                    case '80Hz-Baseline'
                        if ~repeatCycles
                            trialLength=(1/80)*trials.(field).sampleRate*numberCycles;
                            temp=trials.(field).reconstructedPulses_timestamps(pulse_being_processed,1):trialLength:trials.(field).reconstructedPulses_timestamps(pulse_being_processed,2);
                            currentConditionIndex=find(strcmp(conditionIndex,condition));
                            for j=1:length(temp)
                                result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[temp(j) temp(j)+trialLength-1 trials.reconstructedPulses_identities(pulse_being_processed,:)];
                                conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                            end
                        elseif repeatCycles
                            trialLength=(1/80)*trials.(field).sampleRate*numberCycles;
                            temp=trials.(field).reconstructedPulses_timestamps(pulse_being_processed,1):trialLength/numberCycles:trials.(field).reconstructedPulses_timestamps(pulse_being_processed,2);
                            currentConditionIndex=find(strcmp(conditionIndex,condition));
                            for j=1:length(temp)-1
                                result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[temp(j) temp(j)+trialLength-1 trials.reconstructedPulses_identities(pulse_being_processed,:)];
                                conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                            end
                        end
                 end
            end
            pulse_being_processed=pulse_being_processed+1;
        elseif contains(currentCondition,"R") %if current condition is random
            for i=pulse_being_processed:size(trials.reconstructedPulses_identities,1)
                if trials.reconstructedPulses_identities(i,1)~=trials.reconstructedPulses_identities(pulse_being_processed)
                        break;
                end
            end
            temp1=char(regexprep(currentCondition,'-',''));
            for condition={['5.5Hz-' temp1],['40Hz-' temp1],['80Hz-' temp1]}
                switch char(condition)
                    case ['5.5Hz-' char(temp1)]
                        if ~repeatCycles
                            trialLength=(1/5.5)*trials.(field).sampleRate*numberCycles;
                            temp=trials.(field).reconstructedPulses_timestamps(pulse_being_processed,1):trialLength:trials.(field).reconstructedPulses_timestamps(i-1,2);
                            currentConditionIndex=find(strcmp(conditionIndex,condition));
                            for j=1:length(temp)
                                result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[temp(j) temp(j)+trialLength-1 trials.reconstructedPulses_identities(pulse_being_processed,:)];
                                conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                            end
                        elseif repeatCycles
                            trialLength=(1/5.5)*trials.(field).sampleRate*numberCycles;
                            temp=trials.(field).reconstructedPulses_timestamps(pulse_being_processed,1):trialLength/numberCycles:trials.(field).reconstructedPulses_timestamps(i-1,2);
                            currentConditionIndex=find(strcmp(conditionIndex,condition));
                            for j=1:length(temp)-1
                                result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[temp(j) temp(j)+trialLength-1 trials.reconstructedPulses_identities(pulse_being_processed,:)];
                                conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                            end
                        end
                    case ['40Hz-' char(temp1)]
                        if ~repeatCycles
                            trialLength=(1/40)*trials.(field).sampleRate*numberCycles;
                            temp=trials.(field).reconstructedPulses_timestamps(pulse_being_processed,1):trialLength:trials.(field).reconstructedPulses_timestamps(i-1,2);
                            currentConditionIndex=find(strcmp(conditionIndex,condition));
                            for j=1:length(temp)
                                result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[temp(j) temp(j)+trialLength-1 trials.reconstructedPulses_identities(pulse_being_processed,:)];
                                conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                            end
                        elseif repeatCycles
                            trialLength=(1/40)*trials.(field).sampleRate*numberCycles;
                            temp=trials.(field).reconstructedPulses_timestamps(pulse_being_processed,1):trialLength/numberCycles:trials.(field).reconstructedPulses_timestamps(i-1,2);
                            currentConditionIndex=find(strcmp(conditionIndex,condition));
                            for j=1:length(temp)-1
                                result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[temp(j) temp(j)+trialLength-1 trials.reconstructedPulses_identities(pulse_being_processed,:)];
                                conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                            end
                        end
                    case ['80Hz-' char(temp1)]
                        if ~repeatCycles
                            trialLength=(1/80)*trials.(field).sampleRate*numberCycles;
                            temp=trials.(field).reconstructedPulses_timestamps(pulse_being_processed,1):trialLength:trials.(field).reconstructedPulses_timestamps(i-1,2);
                            currentConditionIndex=find(strcmp(conditionIndex,condition));
                            for j=1:length(temp)
                                result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[temp(j) temp(j)+trialLength-1 trials.reconstructedPulses_identities(pulse_being_processed,:)];
                                conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                            end
                        elseif repeatCycles
                            trialLength=(1/80)*trials.(field).sampleRate*numberCycles;
                            temp=trials.(field).reconstructedPulses_timestamps(pulse_being_processed,1):trialLength/numberCycles:trials.(field).reconstructedPulses_timestamps(i-1,2);
                            currentConditionIndex=find(strcmp(conditionIndex,condition));
                            for j=1:length(temp)-1
                                result{currentConditionIndex}(conditionIndexRow{currentConditionIndex},:)=[temp(j) temp(j)+trialLength-1 trials.reconstructedPulses_identities(pulse_being_processed,:)];
                                conditionIndexRow{currentConditionIndex}=conditionIndexRow{currentConditionIndex}+1;
                            end
                        end
                end
            end
            pulse_being_processed=i;
        elseif contains(currentCondition,"occluded") %if current condition is occluded
            break; %ADD CODE FOR OCCLUDED CONDITIONS
        end
    end
end

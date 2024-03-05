%Plots example IED.
%2024/02/26

%% define directories:
root_dir=define_flicker_root_dir;
figures_dir=[root_dir '/stg-analyses/NatureComm2024-figures'];

metadata_tbl=readtable([root_dir '/FlickerStudyMetadata.xlsx'],'Sheet','Subjects','PreserveVariableNames',1);
metadata_tbl=metadata_tbl(:,{'Subject_ID','Subject_ID','Seizure focus broader classification'});

%% fetch data for example session:

examples=readtable([figures_dir '/NatureComm2024_examples.xlsx']);
examples=table2array(examples(strcmp(examples.example,'Figure_6A'),2:end));
example_session=examples(1,1:3); %pick a session to draw example IEDs from

%set some criteria:
sample_rate=200; %sample rate of ied preprocessed data
time_window=0.1; %time window to consider spikes to be part of same event
max_channels=11; %maximum number of channels that can be involved in detecting given event

%define session information:
s=1;
sessions=struct;
sessions(s).subjectID=examples{1};
sessions(s).task=examples{2};
sessions(s).ses=examples{3};

%read ied data and clean:
sessions(s).ied_tbl=readtable([root_dir '/stg-preproc/sub-' sessions(s).subjectID '/task-' sessions(s).task '/ses-' sessions(s).ses '/IED-preproc/sub-' sessions(s).subjectID '_allspikes.csv'],'Delimiter',','); %read detected ieds
sessions(s).ied_tbl(logical(sessions(s).ied_tbl.predicted_class),:)=[]; %remove spikes considered not actual ieds
sessions(s).ied_tbl(:,{'subject','clip_ids','clip','predicted_class'})=[]; %remove info we don't need
sessions(s).ied_tbl.chan=sessions(s).ied_tbl.chan+1; %increment index of channels by 1 (0 in python is 1 in matlab)

%convert spikes to events (i.e. a train of rapidly succeeding spikes, within time_window, is 1 event):
%initialize event table:
sessions(s).event_tbl=table('Size',[0 7],'VariableNames',{'start',...           %start of event, in time samples
    'start_mean',...      %mean start of event (i.e. mean of starts of all spikes constituting the event), in time samples
    'duration',...        %duration of the event (with last spike start + time_window), in time samples
    'spike_starts',...    %start times of each of the spikes in the event, in time samples
    'channels',...        %index of channels implicated in detecting each of the spikes in the event
    'num_spikes',...      %number of spikes in the event
    'num_channels'},...   %number of unique channels implicated in detecting spikes in this event
    'VariableTypes',{'double','double','double','cell','cell','double','double'});
%fill-in event table:
%STOPPED REVIEWING HERE...
spike_list=1; %keeps track of index of spikes to include in given event; start with spike of index 1
for sp=1:size(sessions(s).ied_tbl,1) %for each detected ied
    if sp==size(sessions(s).ied_tbl,1) || sessions(s).ied_tbl.start(sp+1)-sessions(s).ied_tbl.start(sp)>time_window*sample_rate %if we've reached end of ied tbl or the next spike is farther than time window, the current ied event is done and should be added to event table
        sessions(s).event_tbl(end+1,:)={sessions(s).ied_tbl.start(spike_list(1)),...
            mean(sessions(s).ied_tbl.start(spike_list)),...
            ((sessions(s).ied_tbl.start(spike_list(end))+time_window*sample_rate)-sessions(s).ied_tbl.start(spike_list(1))),...
            mat2cell(sessions(s).ied_tbl.start(spike_list)',1,length(spike_list)),...
            mat2cell(sessions(s).ied_tbl.chan(spike_list)',1,length(sessions(s).ied_tbl.chan(spike_list)')),...
            length(spike_list),...
            length(unique(sessions(s).ied_tbl.chan(spike_list)))};
        spike_list=sp+1; %start new list of index of spikes i.e. new event, starting with next spike
    elseif sessions(s).ied_tbl.start(sp+1)-sessions(s).ied_tbl.start(sp)<=time_window*sample_rate %if next spike is within time window
        spike_list=[spike_list,sp+1]; %include the next ied in the current event
    end
end
%sanity check that event table makes sense:

sessions(s).event_tbl(sessions(s).event_tbl.num_channels>max_channels,:)=[]; %remove events that were detected by more than max_channels (considered noise)

%get event start times:
sessions(s).ied_event_times=sessions(s).event_tbl.start(:)'/sample_rate; %convert to absolute times in seconds

%get information about trial times:
sessions(s).trials=importdata([root_dir '/stg-preproc/sub-' sessions(s).subjectID '/task-' sessions(s).task '/ses-' sessions(s).ses '/sub-' sessions(s).subjectID '_stg-preproc_task-' sessions(s).task '_ses-' sessions(s).ses '_nat-beh.mat'],'trials');

%get channel labels for events table:
sessions(s).preproc_labels=readcell([root_dir '/stg-preproc/sub-' sessions(s).subjectID '/task-' sessions(s).task '/ses-' sessions(s).ses '/IED-preproc/sub-' sessions(s).subjectID '_eegdata_labels.csv']);
sessions(s).preproc_labels=sessions(s).preproc_labels';

%ADD CODE TO CHECK THAT ALL PREPROC LABELS MATCH ALL MODULATION
%ANALYSIS LABELS

%fetch subject anat data:
temp=fetch_subject_data(root_dir,{sessions(s).subjectID},'anat');
sessions(s).anat=temp.anat;

%fetch preprocessed LFP data if this is the session we want to show:
sessions(s).data=readtable([root_dir '/stg-preproc/sub-' sessions(s).subjectID '/task-' sessions(s).task '/ses-' sessions(s).ses '/IED-preproc/sub-' sessions(s).subjectID '_eegdata.csv'],'Delimiter',',','PreserveVariableNames',1);
sessions(s).data=table2array(sessions(s).data)';

%% plot IED on LFP:

%plot example detected events, during a stim trial:
s=1;
condition='40Hz-AV'; %pick condition we want to represent
duration=3; %duration of stim trial to present, in seconds
flicker_color=condition_color(condition); %determine modality color

trial_nber=115; %pick one of those trials
trial_timestamps=sessions(s).trials.clinrecording.trials_timestamps(trial_nber,:); %get start and end samples for that trial
trial_timestamps=trial_timestamps/sessions(s).trials.clinrecording.sampleRate; %convert timestamps to seconds
trial_timestamps(2)=trial_timestamps(1)+3;

%determine list of events that occurred during picked time period:
event_list=sessions(s).event_tbl(sessions(s).event_tbl.start>=trial_timestamps(1)*sample_rate & sessions(s).event_tbl.start<=trial_timestamps(2)*sample_rate,:); %find events that occurred during that trial
depths_of_interest=unique([event_list.channels{:}]);
depths_of_interest=sessions(s).preproc_labels(depths_of_interest);
depths_of_interest=extract_clinLFP_labels(depths_of_interest);
depths_of_interest={depths_of_interest.depth_electrode_name};
contacts_of_interest=arrayfun(@(x) sessions(s).preproc_labels(startsWith(sessions(s).preproc_labels,x)),depths_of_interest,'UniformOutput',false);
contacts_of_interest=[contacts_of_interest{:}]'; %list of contacts we will need to plot

nber_timeseries_perplot=length(contacts_of_interest); %number of timeseries we're going to plot

%create set of colors to describe different events:
if time_window<0
    event_colors=repmat([1 0 0],size(event_list,1),1);
else
    event_colors=rand([size(event_list,1),3]);
end

%plot timeseries for each channel:
figure_making('height',4,'width',7);
hold on;
space_btw_timeseries=0.0007; %400microV (?)
ylim([0-(nber_timeseries_perplot-1)*space_btw_timeseries-space_btw_timeseries/2-length(depths_of_interest)*space_btw_timeseries-space_btw_timeseries 0+space_btw_timeseries/2+space_btw_timeseries]);
temp=gca;
temp.YTickLabel=[];
temp.YTick=[];
temp.YAxis.Visible='off';
temp_ylim=temp.YLim;
temp_ylim(2)=temp_ylim(2)+range(temp_ylim)/100;
add_scale_bars(gca,0,0.001,0.01);
plot_flicker(str2num(regexprep(condition,'Hz.+','')),temp.YLim,duration,flicker_color,0,100);
ylim(temp_ylim);
row_nber=0;
depth_electrode_name_added=0;
for ch=1:length(contacts_of_interest) %for each channel
    
    channel_nber=find(strcmp(sessions(s).preproc_labels,contacts_of_interest(ch)));
    
    if ch>1
        set(gca, 'XLimSpec', 'Tight');
    end
    
    plot(0:1/sample_rate:duration,sessions(s).data(channel_nber,round(trial_timestamps(1)*sample_rate):round(trial_timestamps(2)*sample_rate))-space_btw_timeseries*row_nber,'k');
    temp1=extract_clinLFP_labels(contacts_of_interest(ch));
    temp2=regexp(temp1.channel_names,'\d+');
    temp2=temp1.channel_names{1}(temp2{:}(end):end);
    %text(-0.15,mean(data(channel_nber,round(trial_timestamps(1)*sample_rate):round(trial_timestamps(2)*sample_rate)))-space_btw_timeseries*row_nber,temp2,'FontSize',6);
    
    if ~depth_electrode_name_added
        num_contacts=sum(startsWith(contacts_of_interest,temp1.depth_electrode_name));
        text(-duration/15,-space_btw_timeseries*row_nber-(num_contacts*space_btw_timeseries)/2,temp1.depth_electrode_name);
        depth_electrode_name_added=1;
    end
    
    %highlight any spikes that occured in that channel, with colors
    %corresponding to individual events:
    temp=find(arrayfun(@(x) any(ismember(x{:},channel_nber)),event_list.channels)); %find events in list that were picked up by that channel
    if ~isempty(temp) %if some of events in list were picked up by that channel
        for i=temp' %for each of those picked up events
            line_color=event_colors(i,:); %determine color of event
            for j=find(ismember(event_list.channels{i},channel_nber)) %for each channel that picked up that event
                x=(event_list.spike_starts{i}(j)-round(trial_timestamps(1)*sample_rate))/sample_rate:1/sample_rate:(event_list.spike_starts{i}(j)-round(trial_timestamps(1)*sample_rate))/sample_rate+0.1; %define x values for this spike
                y=sessions(s).data(channel_nber,event_list.spike_starts{i}(j):event_list.spike_starts{i}(j)+sample_rate*0.1)-space_btw_timeseries*row_nber;
                y_patch=[repmat(min(y),1,length(y)) repmat(max(y),1,length(y))];
                %p=patch([x x(end:-1:1)],y_patch,line_color,'EdgeColor','none','FaceAlpha',0.6);
                %uistack(p,'bottom');
                plot(x,y,'Color','r','LineWidth',1);
            end
        end
    end
    
    if ~(ch+1>=length(contacts_of_interest))
        temp_depth_electrode=extract_clinLFP_labels(contacts_of_interest(ch+1));
        temp_depth_electrode=temp_depth_electrode.depth_electrode_name;
        if ~startsWith(contacts_of_interest(ch),temp_depth_electrode)
            row_nber=row_nber+2;
            depth_electrode_name_added=0;
        else
            row_nber=row_nber+1;
        end
    else
        row_nber=row_nber+1;
    end
end

xlabel('Time (s)');
ylabel(['LFP contact' newline newline newline]);

figure_making('operation','save','filename',[figures_dir '/Figure_6A.pdf']);

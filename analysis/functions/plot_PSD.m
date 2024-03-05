%Plots PSD for given flicker condition and channel.
%2024/02/25

function plot_PSD(condition,channel,PSD_results,channel_labels,condition_code,color,std_error,show_trials,hide_60Hz)
    %plot PSD:
    psd_result=PSD_results.data{strcmp(channel_labels,channel),strcmp(condition_code,condition)}{1};
    plot(PSD_results.data{strcmp(channel_labels,channel),strcmp(condition_code,condition)}{3},log10(mean(psd_result)),'Color',color);
    
    %add standard error and traces for individual trials if requested:
    if std_error
        x=[PSD_results.data{strcmp(channel_labels,channel),strcmp(condition_code,condition)}{3} PSD_results.data{strcmp(channel_labels,channel),strcmp(condition_code,condition)}{3}(end:-1:1)];
        p=patch(x,log10([mean(psd_result)-std(psd_result)/sqrt(size(psd_result,1)) mean(psd_result(:,end:-1:1))+std(psd_result(:,end:-1:1))/sqrt(size(psd_result,1))]),color,'FaceAlpha',0.2,'EdgeColor','none');
        %p=patch(x,log10([mean(psd_result)-std(psd_result) mean(psd_result(:,end:-1:1))+std(psd_result(:,end:-1:1))]),color,'FaceAlpha',0.2,'EdgeColor','none');
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    elseif show_trials
        for i=1:size(psd_result,1)
            hold on;
            p=plot(PSD_results.data{strcmp(channel_labels,channel),strcmp(condition_code,condition)}{3},log10(psd_result(i,:)),'Color',color,'LineWidth',1);
            p.Color(4)=0.1;
        end
    end
    
    %hide ground noise signal:
    if hide_60Hz
        p=gca;
        min_val=min(log10(mean(PSD_results.data{strcmp(channel_labels,channel),strcmp(condition_code,condition)}{1})));
        max_val=max(log10(mean(PSD_results.data{strcmp(channel_labels,channel),strcmp(condition_code,condition)}{1})));
        p=patch([58 62 62 58],[min_val min_val max_val max_val],'w','EdgeColor','none','FaceAlpha',1);
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end 
end

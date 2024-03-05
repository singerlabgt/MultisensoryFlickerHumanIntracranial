%Plots PSTH showing response of given neuronal unit to senosry flicker for
%given condition.
%2024/02/25

function plot_flicker_psth(psth_data,conditionIndex,condition,unit_label,max_spikerate)
    %determine condition color:
    temp=strsplit(condition,'-');
    if strcmp(temp{2},'V') %visual in orange
        condition_color=[0.9100 0.4100 0.1700];
    elseif strcmp(temp{2},'A') %auditory in blue
        condition_color='b';
    elseif strcmp(temp{2},'AV') %audio-visual in green
        condition_color='g';
    end
    
    if any(psth_data{strcmp(conditionIndex,condition)}.avg(strcmp(psth_data{strcmp(conditionIndex,condition)}.label,unit_label),:)) %only plot PSTH if we have any spike rate for that unit and condition
        %plot PSTH for that stim condition using FieldTrip:
        cfg=[];
        cfg.errorbars='no';
        cfg.spikechannel=unit_label;
        ft_spike_plot_psth(cfg,psth_data{strcmp(conditionIndex,condition)});
        hold on;
        %clean and color plot:
        h = get(gca, 'Children');
        set(h,'FaceColor',condition_color)
        set(gca,'XTickLabel',[]);
        set(gca,'xtick',[]);
        set(gca,'XLabel',[]);
        set(gca,'YLabel',[]);

        %plot PSTH for the comparator condition using FieldTrip:
        temp=strsplit(condition,'-');
        comparator=[temp{1} '-R' temp{2}];
        ft_spike_plot_psth(cfg,psth_data{strcmp(conditionIndex,comparator)});
        %clean and color plot:
        h = get(gca, 'Children');
        set(h(1), 'FaceColor', [0.5 0.5 0.5]);
        h(1).YData=h(1).YData*(-1);
        set(gca,'YLabel',[]);
        set(gca,'XLabel',[]);

        %add labels of interest:
        freq=str2double(regexprep(temp{1},'Hz',''));
        temp=gca;
        temp_xticks=0:1/freq/2:temp.XLim(2);
        xticks(temp_xticks);
        xticklabels(cellstr(strtrim(string(num2str([temp_xticks*1000]'))))');
        xlabel('Time (ms)');

        %set Y-lims:
        if isnan(max_spikerate)
            Xlim=xlim;
            axis tight;
            Ylim=ylim;
            ylim([Ylim(1)-(Ylim(2)-Ylim(1))*0.05 Ylim(2)+ (Ylim(2)-Ylim(1))*0.05]);
            xlim(Xlim);
        elseif ~isnan(max_spikerate)
            maxYaxisValue=max_spikerate;
            set(gca,'ylim',[-maxYaxisValue-maxYaxisValue*0.1 maxYaxisValue+maxYaxisValue*0.1]);
        end
        temp=xlim;
        xlim([0 ceil(temp(2)*1000)/1000]);
        temp=ylim;
        plot_flicker(freq,temp,0.5,condition_color,0);
        temp=gca;
        ylim([temp.YLim(1)-((temp.YLim(2)-temp.YLim(1))/20) temp.YLim(2)+((temp.YLim(2)-temp.YLim(1))/20)]);
        temp=get(gca,'YTickLabel');
        set(gca,'YTickLabel',cellfun(@(x) num2str(abs(str2double(x))),temp,'UniformOutput',false));
    end
    
    ylabel('Spike rate');
    title(unit_label,'Interpreter','none');
end

%assigns color based on flicker modulation amplitude.
%2024/02/25

function colors=assign_zscore_color(zscores,threshold)
    %create colors from yellow (low modulation) to red (high modulation):
    plot_colormap=[flipud(autumn)];
    
    %assign colors:
    colors=[];
    for i=1:length(zscores) %for each modulation value
        if zscores(i)<=0 %if modulation value is negative or 0, assign grey color
            colors(i,:)=[0.5 0.5 0.5];
        elseif zscores(i)>0 && zscores(i)<=threshold %if modulation value is in-between 0 and capped max value, assign relative yellow-to-red color
            colors(i,:)=vals2colormap(zscores(i),plot_colormap,[0 threshold]);
        elseif zscores(i)>threshold %if modulation value is greater than capped max value, assign max red value
            colors(i,:)=plot_colormap(end,:);
        end
    end
end

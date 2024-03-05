

function plot_ERP(ERP_results,condition,channel,std_error,show_trials,show_stim)

    %find which ERP_results correspond to condition
    cond=find(strcmp(ERP_results{end},condition));

    if endsWith(condition,'occluded_AV') || strcmp(condition,'Baseline')
        line_color='k';
    elseif endsWith(condition,'AV')
        line_color='g';
    elseif endsWith(condition,'V')
        line_color=[0.9100 0.4100 0.1700];
    elseif endsWith(condition,'A')
        line_color='b';
    end
    
    %get erp results:
    erp_result=reshape(ERP_results{cond}.trial(:,strcmp(ERP_results{cond}.label,channel),:),size(ERP_results{cond}.trial,1),size(ERP_results{cond}.trial,3));
    
    %plot erp results:
    plot(ERP_results{cond}.time,mean(erp_result),'Color',line_color,'LineWidth',0.5);
              
    %plot standard error if requested:
    if std_error
        x=[ERP_results{cond}.time ERP_results{cond}.time(end:-1:1)];
        p=patch(x,[mean(erp_result)-std(erp_result)/sqrt(size(erp_result,1)) mean(erp_result(:,end:-1:1))+std(erp_result(:,end:-1:1))/sqrt(size(erp_result,1))],line_color,'FaceAlpha',0.1,'EdgeColor','none');
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    
    %show trials if requested:
    if show_trials
        
    end
    
    yl=ylim;
    axis tight;
    ylim(yl);
    
    if show_stim
        temp1=get(gca,'ylim');

        x = [0; 0.0125; 0.0125; 0];
        y = [temp1(1); temp1(1); temp1(2); temp1(2)];

        p=patch(x,y,[0.5 0.5 0.5],'EdgeColor','none','FaceAlpha','0.5');
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        xline(0,'-','','LabelHorizontalAlignment','left');
    end
    
    hold on;
    
    xlabel('Time (s)');
    ylabel('Electric Potential (uV)');
    set(get(gca,'YLabel'),'Rotation',90);
    plot([ERP_results{cond}.time(1), ERP_results{cond}.time(end)], [0 0], 'k--'); % add horizontal line
    set(gca,'children',flipud(get(gca,'children')));
end




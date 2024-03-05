%Plots PSDs for all flicker conditions (except random and occluded),
%all electrodes (1 depth electrode per figure), and saves the figure. 
%The following colors correspond to:
%*orange: visual modality.
%*blue: auditory modality.
%*green: audio-visual modality.
%2024/02/25

function plot_depthelectrode_PSDs(PSD_results,outputFolder)
    depth_electrodes=extract_clinLFP_labels(PSD_results.label); %organize electrode contacts by detph electrode
    
    fig=[];
    for i=1:length(depth_electrodes) %for each detph electrode
        fig(i)=figure('units','normalized','outerposition',[0 0 1 1]);
        tiledlayout(4,5);
        for j=1:length(depth_electrodes(i).channel_names) %for each electrode contact
            nexttile;
            title(depth_electrodes(i).channel_names{j});
            hold on;
            for condition=1:length(PSD_results.condition) %for each condition
                if ~isempty(PSD_results.condition{condition}) && ~contains(PSD_results.condition{condition},'occluded')
                    %determine condition color:
                    temp=strsplit(PSD_results.condition{condition},'-');
                    if length(temp)==1 %baseline condition in black
                        line_color='k';
                    else
                        switch temp{2} %stimulation conditions in their respective colors
                            case 'V'
                                line_color=[0.9100 0.4100 0.1700]; %visual in orange
                            case 'A'
                                line_color='b'; %auditory in blue
                            case 'AV'
                                line_color='g'; %audio-visual in green
                        end
                    end
                    
                    %plot PSD result:
                    psd_result=PSD_results.data{strcmp(PSD_results.label,depth_electrodes(i).channel_names{j}),condition}; %get PSD results of interest
                    plot(psd_result{3},log10(mean(psd_result{1})),'Color',line_color,'LineWidth',1); %plot
                    
%                     %plot standard error:
%                     x=[psd_result{3} psd_result{3}(end:-1:1)];
%                     p=patch(x,log10([mean(psd_result{1})-std(psd_result{1})/sqrt(size(psd_result{1},1)) mean(psd_result{1}(:,end:-1:1))+std(psd_result{1}(:,end:-1:1))/sqrt(size(psd_result{1},1))]),line_color,'FaceAlpha',0.2,'EdgeColor','none');
%                     set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');       
                end
            end
            
            %clean plot:
            p=gca;
            p=patch([58 62 62 58],[p.YLim(1) p.YLim(1) p.YLim(2) p.YLim(2)],'w','EdgeColor','none','FaceAlpha',0.9); %hide ground noise peak
            set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            xline(5.5,'--');
            xline(40,'--');
            xline(80,'--');
        end
        
        %print figure:
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,...
            'PaperPosition',[0 0 screenposition(3:4)],...
            'PaperSize',[screenposition(3:4)]);
        print(fig(i),[outputFolder '/depth-electrode-' depth_electrodes(i).depth_electrode_name '_entrainment-PSD.pdf'],'-dpdf','-fillpage');
    end
    close all;
end

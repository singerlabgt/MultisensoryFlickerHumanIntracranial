%Plots flicker stimulus on plot. 
%2024/02/25

function plot_flicker(freq,y_lim,duration,patch_color,single_pulse,varargin)
    %COMMENTS EXPLAINING BELOW?
    if ~isempty(varargin)
        height_divider=varargin{1};
    else
        height_divider=5;
    end
    yrange=[y_lim(2) y_lim(2)+((y_lim(2)-y_lim(1))/height_divider)];
    start_end=[0 1/freq/2;yrange];
    
    if single_pulse %if ploting single pulse (for spep task)
        %plot single pulse in condition color:
        p=patch([start_end(1,:) start_end(1,end:-1:1)],[start_end(2,1) start_end(2,1) start_end(2,2) start_end(2,2)],patch_color,'EdgeColor','none');
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        start_end(1,:)=start_end(1,:)+1/freq/2;
        
        %plot rest (no stim) in black
        p=patch([start_end(1,1) duration duration start_end(1,1)],[start_end(2,1) start_end(2,1) start_end(2,2) start_end(2,2)],'k','EdgeColor','none');
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
    elseif ~single_pulse %if not plotting single pulse (for flicker task)
        for i=1:(duration/(1/freq)) %for duration requested
            %plot pulses in condition color
            p=patch([start_end(1,:) start_end(1,end:-1:1)],[start_end(2,1) start_end(2,1) start_end(2,2) start_end(2,2)],patch_color,'EdgeColor','none');
            set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            start_end(1,:)=start_end(1,:)+1/freq/2;

            %plot space between pulses in black
            p=patch([start_end(1,:) start_end(1,end:-1:1)],[start_end(2,1) start_end(2,1) start_end(2,2) start_end(2,2)],'k','EdgeColor','none');
            set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            start_end(1,:)=start_end(1,:)+1/freq/2;
        end
    end
end

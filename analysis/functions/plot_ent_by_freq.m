%Plots modulation or PLV amplitude (and significance) by flicker frequency of stimulation.
%2024/02/26

function plot_ent_by_freq(zscore_table,pvalue_table,ch,display_pvalue,varargin)
    if length(varargin)==1 && strcmp(varargin{1},'plv')
        line_style=':';
        marker_style='d';
    else
        line_style='-';
        marker_style='o';
    end

    temp=arrayfun(@(x) string(strsplit(x{:},'-')),zscore_table.Properties.VariableNames','UniformOutput',false);

    %determine modality we're analyzing, and pick corresponding plotting color:
    modality=temp{1,1}{2};
    line_color=condition_color(modality);
    
    %plot amplitudes:
    frequencies=arrayfun(@(x) str2num(regexprep(x{1}{1,1},'Hz','')),temp,'UniformOutput',false);
    frequencies=[frequencies{:}];
    [~,freq_index]=sort(frequencies,'ascend'); %sort stim frequencies in ascending order
    plot(frequencies(freq_index),zscore_table{ch,freq_index},'Color',line_color,'LineStyle',line_style);
    hold on;
    
    %plot significance:
    for condition=zscore_table.Properties.VariableNames %for each condition
        frequency=strsplit(condition{:},'-');
        frequency=str2num(regexprep(frequency{1},'Hz',''));
        if ~display_pvalue
            scatter(frequency,zscore_table{ch,condition},'filled','MarkerFaceAlpha',line_color,'Marker',marker_style);
        elseif display_pvalue
            if pvalue_table{ch,condition}<0.05
                alpha_value=1; %significant result is opaque
            else
                alpha_value=0.2; %non-significant result is transparent
            end
            scatter(frequency,zscore_table{ch,condition},'filled','MarkerFaceColor',line_color,'MarkerFaceAlpha',alpha_value,'Marker',marker_style);
        end
        
        %set X-lim:
        if strcmp(condition,zscore_table.Properties.VariableNames{1})
            hold on;
            xlim([0 85.5]);
        end
    end
    
    %set plot labels:
    title(ch);
    xlabel('Stimulation frequency');
    ylabel('Fold change in power');
end

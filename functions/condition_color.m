%Returns the color the condition should be (visual is orange, audio is
%blue, audio-visual is green).
%2024/02/25

function line_color=condition_color(condition)
    if iscell(condition)
        condition=condition{:};
    end
    condition=regexprep(condition,'.+-','');
    switch condition
        case 'V'
            line_color=[0.9100 0.4100 0.1700]; %orange
        case 'A'
            line_color=[0 0 1]; %blue
        case 'AV'
            line_color=[0 1 0]; %green
        case 'occluded_AV'
            line_color=[0 0 0]; %black
    end
end

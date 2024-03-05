%Orders flicker conditions in increasing frequency order.

function conditions=order_flicker_conditions(conditions)
    temp=arrayfun(@(x) str2double(regexprep(x,'Hz.+','')),conditions);
    [~,temp]=sort(temp);
    conditions=conditions(temp);
end
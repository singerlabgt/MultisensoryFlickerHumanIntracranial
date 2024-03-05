

function asterisk=pval_to_asterisk(pval)
    if pval<0.001
        asterisk='***';
    elseif pval<0.01
        asterisk='**';
    elseif pval<0.05
        asterisk='*';
    else
        asterisk='';
    end
end
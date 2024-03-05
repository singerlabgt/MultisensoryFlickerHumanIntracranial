
function add_scale_bars(ax,quantity_xaxis,quantity_yaxis,shift_from_origin)
    %x scale bar:
    disp(['x scale bar is ' num2str(quantity_xaxis)]);
    annotation('line',[ax.Position(1) ax.Position(1)+(quantity_xaxis*ax.Position(3))/range(ax.XLim)]-shift_from_origin,[ax.Position(2) ax.Position(2)]-shift_from_origin,'Linewidth',1);
    
    %y scale bar:
    disp(['y scale bar is ' num2str(quantity_yaxis)]);
    annotation('line',[ax.Position(1) ax.Position(1)]-shift_from_origin,[ax.Position(2) ax.Position(2)+(quantity_yaxis*ax.Position(4))/range(ax.YLim)]-shift_from_origin,'Linewidth',1);
end
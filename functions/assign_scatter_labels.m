%Assigns labels provided by labels_array to each dot in the scatter plot.

function assign_scatter_labels(s,labels_array)
    row = dataTipTextRow('Label',labels_array);
    s.DataTipTemplate.DataTipRows(end+1) = row;
end
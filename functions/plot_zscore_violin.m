

function v=plot_zscore_violin(vectors_array,names_array,threshold)
    
    num_groups=length(vectors_array);
    vectors=[];
    names={};
    for i=1:num_groups
        vectors=[vectors;vectors_array{i}];
        names=[names; repmat(names_array(i),[1 size(vectors_array{i},1)])'];
    end
    
    %note: function below is from Violinplot-Matlab (not FieldTrip)
    v=violinplot(vectors,names,'EdgeColor',[0 0 0],'GroupOrder',names_array);
    for i=1:num_groups
        v(i).ViolinAlpha=0;
        v(i).ViolinPlot.EdgeAlpha=0.3;
        v(i).ScatterPlot.MarkerFaceAlpha=0.5;
        %colors=assign_zscore_color(vectors_array{i},threshold);
        v(i).ScatterPlot.MarkerFaceColor='flat';
        v(i).ScatterPlot.CData=repmat([0 0 0],[length(vectors_array{i}),1]);
        v(i).ScatterPlot.SizeData=5;
        v(i).BoxPlot.FaceAlpha=0;
        v(i).BoxPlot.EdgeAlpha=0;
        %v(i).WhiskerPlot.Visible='off';
    end
    
end
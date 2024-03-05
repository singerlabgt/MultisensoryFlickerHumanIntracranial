%this function is used to create and save figures with set formats for
%paper, poster or PPT figures:

function figure_making(NameValueArgs)

    arguments
        NameValueArgs.operation char {mustBeMember(NameValueArgs.operation,{'create','save'})} = 'create'
        NameValueArgs.format char {mustBeMember(NameValueArgs.format,{'paper','poster','ppt'})} = 'paper'
        NameValueArgs.width double=1
        NameValueArgs.height double=1.2
        NameValueArgs.filename char
    end
    
    
    
    if strcmp(NameValueArgs.operation,'create') %means we want to create figure
        
        %set font:
        set(groot,'defaultAxesFontName','Arial')
        set(groot,'defaultAxesFontSizeMode','manual');
        
        %determine by how much font sizes should be multiplied:
        if strcmp(NameValueArgs.format,'paper') %font multiplier for figure
            multiplier=1;
        elseif strcmp(NameValueArgs.format,'poster') %font multiplier for poster
            multiplier=2;
        elseif strcmp(NameValueArgs.format,'ppt') %font multiplier for ppt
            multiplier=3;
        end
        
        %set font sizes:
        set(groot,'defaultAxesFontSize',6*multiplier);
        set(groot,'defaultAxesLabelFontSizeMultiplier',7/6);
        set(groot,'defaultTextFontSize',6*multiplier*7/6);
        set(groot,'defaultTextarrowshapeFontSize',6*multiplier*7/6);
        set(groot,'defaultTextboxshapeFontSize',6*multiplier*7/6);
        set(groot,'defaultColorbarFontSize',6*multiplier*7/6);
        set(groot,'defaultAxesTitleFontSizeMultiplier',12/6);
        set(groot,'defaultAxesLineWidth',0.5*multiplier);
        set(groot,'defaultLineLineWidth',0.5*multiplier);
        set(groot,'defaultLineMarkerSize',3*multiplier);
        
        
        %create figure of given height and width:
        figure('Units','inches','Position',[0.5 0.5 NameValueArgs.width  NameValueArgs.height],'PaperUnits','inches','PaperSize',[NameValueArgs.width NameValueArgs.height],'PaperPosition',[0 0 NameValueArgs.width NameValueArgs.height]); %make figure with proper dimensions
    
    elseif strcmp(NameValueArgs.operation,'save')
        
        h=findobj(gcf, '-property','ZData');
        h=get(h,'ZData');
        if any(~arrayfun(@(x) isempty(x{:}),h))
            disp('Found 3D plot- using export_fig function...');
            set(gcf,'Color','w');
            export_fig(NameValueArgs.filename);
        else
            set(gcf,'renderer','painters');
            saveas(gcf,NameValueArgs.filename,'pdf');
        end
        close(gcf);
        
    end
end
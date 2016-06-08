function [] = save_figure(figName, handle, typeFigure, res)
% function [] = save_figure(figName, handle, typeFigure, res)
% function [] = save_figure(figName, [], [], [])
% Save the current figure
%
%
    if isempty(handle)
        handle = gcf;
    end
    
    figure(handle);
    
    if isempty(res)
        res = 600;
    end
    if isempty(typeFigure)
        typeFigure = 'png';
    end
    
    switch  typeFigure
        %case 'pdf'
           %print([figName '.pdf'], res, '-dpdf', handle);
        case 'png'
           %print([figName '.png'], res, '-dpng', handle);
           print([figName '.png'], '-dpng', '-r300');
        otherwise
           error('Figure type is  png.')
    end
    
end

%    figure(handle)
%    handle = gcf
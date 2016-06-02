function ynew = fakeQuery(yold, ratioStartEnd, shorten)


    if exist('shorten', 'var')
        iStart = round(numel(yold)*ratioStartEnd(1));
        iEnd   = numel(yold)- round(numel(yold)*ratioStartEnd(2));
        y = yold(iStart:iEnd);        
    else
        yStart = yold(1).*ones(1, round(numel(yold)*ratioStartEnd(1)  ));
        yEnd   = yold(end).*ones(1, round(numel(yold)*ratioStartEnd(2)));
        y = [ yStart  yold yEnd];
    end
    
    ynew = interp1(linspace(0,1,numel(y)), y,  linspace(0,1,numel(yold))  );
%     
%     figurew
%     plot( ynew)
%     plot( yold, 'r')
    
    

end
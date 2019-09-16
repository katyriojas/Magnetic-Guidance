function [hp,Xthresh,Yperc,Xinset,Yinset] = plotCDF(X,color,ls,threshold,upperThresh,lowerThresh)

    Xplot = sort(X);
    Yplot = (1:length(X)).'/length(X);
    Yplot = [0; Yplot(1:(end-1))];
    hp = plot(Xplot, Yplot,'Color',color,'LineStyle',ls,'LineWidth',1);
    indx = find(Xplot>threshold,1);
    Xthresh = Xplot(indx);
    Yperc = Yplot(indx);
    
    idxThresh_L = find(Xplot>lowerThresh,1);
    idxThresh_U = find(Xplot>upperThresh,1);
    Xinset = Xplot(idxThresh_L:idxThresh_U);
    Yinset = Yplot(idxThresh_L:idxThresh_U);

end
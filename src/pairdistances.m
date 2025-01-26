function sd = pairdistances(tmpX_all)
 %calculate pairwise distances
    sd=0*ones(size(tmpX_all,2),size(tmpX_all,2));
    for iX=1:size(tmpX_all,2)
        for iY=iX+1:size(tmpX_all,2)
            iarr = ~isnan(tmpX_all(:, iX)) & ~isnan(tmpX_all(:, iY));
            sd(iX,iY)=mean(power(tmpX_all(iarr,iX)-tmpX_all(iarr,iY),2));
            sd(iY,iX)=sd(iX,iY);
        end
    end
end


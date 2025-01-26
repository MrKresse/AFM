function sd2 = logdistances(tmpX_all,sd)
    sd2=0*ones(size(tmpX_all,2),size(tmpX_all,2));
    for iX=1:size(tmpX_all,2)
        for iY=iX+1:size(tmpX_all,2)
            iarr=log(sd(:,iX))-log(sd(:,iY));
            iarr=iarr(isinf(iarr)==0);
            sd2(iX,iY)=mean(iarr.*iarr);
            sd2(iY,iX)=sd2(iX,iY);
        end
    end
end


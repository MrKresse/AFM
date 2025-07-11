function ii2 = greedycluster(sd2,ii,threshold)
    ii_tmp=ii;
    ii2=[];
    %greedy clustering
    while ~isempty(ii)
        %remove last curve
        ii2=[ii(end) ii2];
        ii=ii(1:end-1);
        iarr1=ii_tmp==ii2(1);
        %find curves closer than the threshold to last curve
        iarr2=find(sd2(iarr1,:)<threshold);
    
        for iX=1:length(iarr2)
            %add similar curves to cluster
           iarr3=find(ii==ii_tmp(iarr2(iX)));
           if ~isempty(iarr3)
               ii2=[ii(iarr3) ii2];
               ii(iarr3)=NaN;
               ii=ii(~isnan(ii));
           end
        end
    end
end


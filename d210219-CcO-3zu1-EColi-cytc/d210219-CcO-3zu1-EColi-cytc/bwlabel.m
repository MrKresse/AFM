function [pixel_map,cent_pos]=bwlabel(bin_map)

    pixel_map=0*bin_map;
    cent_pos=[];

    [cc,hc]=contour(bin_map,1);
    close gcf;
    
    iC=1;
    iP=1;
    while iC<length(cc)
        cont_X=cc(2,iC+1:iC+cc(2,iC));
        cont_Y=cc(1,iC+1:iC+cc(2,iC));
        centr_X=[];
        centr_Y=[];
        
        Xmin=min(round([cont_X]))-1;
        Xmax=max(round([cont_X]))+1;
        
        for iX=Xmin:Xmax
            iarr=find(round(cont_X)==iX);
            if length(iarr)>0
                Ymin=min(round([cont_Y(iarr)]));
                Ymax=max(round([cont_Y(iarr)]));
                pixel_map(iX,Ymin:Ymax)=iP*bin_map(iX,Ymin:Ymax);
                centr_X=[centr_X iX*bin_map(iX,Ymin:Ymax)];
                centr_Y=[centr_Y (Ymin:Ymax).*bin_map(iX,Ymin:Ymax)];
            end
        end
        
        centr_X=centr_X(centr_X>0);
        centr_Y=centr_Y(centr_Y>0);
        cent_pos=[cent_pos; mean(centr_X) mean(centr_Y)];
        
        iC=iC+1+cc(2,iC);
        iP=iP+1;
    end  
end

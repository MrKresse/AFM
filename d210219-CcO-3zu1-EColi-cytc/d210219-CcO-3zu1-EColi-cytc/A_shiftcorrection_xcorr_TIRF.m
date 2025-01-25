clear;

base_str='210219-CcO-3zu1-EColipolar-cytc_0uM-001.nd2.1.tif';

tmp=imreadBFmeta(strcat('./incoming/',base_str));
Ndim=tmp.nframes;
xdim=tmp.width;
ydim=tmp.height;

delta_x=0;
delta_y=0;
delta_max=15;
Nborder=100;

%align 405
tmp_pre=double(imread(strcat('./incoming/',base_str),1));
tmp_pre=reshape(double(tmp_pre)',xdim,ydim);
tmp_pre=interp2(double(tmp_pre),1);
tmp=0*ones(size(tmp_pre,1)+2*Nborder,size(tmp_pre,2)+2*Nborder);
tmp((Nborder+1:Nborder+size(tmp_pre,1))+delta_x,(Nborder+1:Nborder+size(tmp_pre,2))+delta_y)=tmp_pre;
imwrite(uint16(tmp'),strcat('./incoming/',base_str(1:end-4),'.nf.tif'),'WriteMode','overwrite');

for iX=2:Ndim
% for iX=2:5
    display(strcat('filtering405: ',num2str(iX),'; d_x=',num2str(delta_x),'; d_y=',num2str(delta_y)));
    
    tmp_act=double(imread(strcat('./incoming/',base_str),iX));
    tmp_act=reshape(double(tmp_act)',xdim,ydim);
    tmp_act=interp2(double(tmp_act),1);
    
    tmp_corr_mat=0*ones(2*delta_max+1);
    
    CX_arr=(-delta_max:delta_max)+delta_x;
    CY_arr=(-delta_max:delta_max)+delta_y;
    
    for iCX=1:length(CX_arr)
        for iCY=1:length(CY_arr)
            tmp_corr=tmp_pre(200:end-200,200:end-200).*tmp_act(200+CX_arr(iCX):end-200+CX_arr(iCX),200+CY_arr(iCY):end-200+CY_arr(iCY));
            tmp_corr_mat(iCX,iCY)=sum(sum(tmp_corr));
        end
    end
    
    tmp_xcorr=tmp_corr_mat;
    
    [yyx,yix]=max(tmp_xcorr,[],1);
    [yyy,yiy]=max(yyx);
    delta_y=CY_arr(yiy);
    delta_x=CX_arr(yix(yiy));
    
    tmp=0*ones(size(tmp_act,1)+2*Nborder,size(tmp_act,2)+2*Nborder);
    tmp((Nborder+1:Nborder+size(tmp_act,1))-delta_x,(Nborder+1:Nborder+size(tmp_act,2))-delta_y)=tmp_act;
    
    if iX==1
        imwrite(uint16(tmp'),strcat('./incoming/',base_str(1:end-4),'.nf.tif'),'WriteMode','overwrite');
    else
        imwrite(uint16(tmp'),strcat('./incoming/',base_str(1:end-4),'.nf.tif'),'WriteMode','append');
    end
end
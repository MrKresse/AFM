clear;

base_str='210219-CcO-3zu1-EColipolar-cytc_0_05uM-002.nd2.1.nf-blur.tif';
int_thres=126;
% base_str='210219-CcO-3zu1-EColipolar-cytc_0_5uM-004.nd2.1.nf-blur.tif';
% int_thres=126;
% base_str='210219-CcO-3zu1-EColipolar-cytc_0_15uM-003.nd2.1.nf-blur.tif';
% int_thres=126;

pp=[0.001 -0.027 0.247 -1.038 8.204];

tmp=imreadBFmeta(strcat('./incoming/',base_str));
xdim=tmp.width;
ydim=tmp.height;
Ndim=tmp.zsize;
if Ndim>1
    use_zsize=1;
else
    Ndim=tmp.nframes;
    use_zsize=0;
end

if use_zsize==1
    HCO_stack=imreadBF(strcat('./incoming/',base_str),1:Ndim,1,1);
else
    HCO_stack=imreadBF(strcat('./incoming/',base_str),1,1:Ndim,1);
end

HCO_stack_m=median(HCO_stack,3);
BW=(HCO_stack_m>int_thres);
[ves_pix,ves_pos]=bwlabel(BW);

I_back_arr=[];
BW=ves_pix;

for iN=1:Ndim
    HCO_stack_N=HCO_stack(:,:,iN);
    BW_tmp=(HCO_stack_N>int_thres);
    [ves_pix,ves_pos_nil]=bwlabel(BW_tmp);
    BW(:,:,iN)=ves_pix;
    BW_back=(BW_tmp==0);
    I_back=HCO_stack_N(BW_back);
    I_back=I_back(I_back>0);
    I_back_arr=[I_back_arr median(I_back(isnan(I_back)==0))];
end

part_int=[];

for iP=1:size(ves_pos,1)
    iX=floor(ves_pos(iP,1));
    iY=floor(ves_pos(iP,2));
    
    part_int_tmp=[];
    
    for iN=1:Ndim
        HCO_stack_N=HCO_stack(:,:,iN);
        BW_tmp=BW(:,:,iN);
        part_ID=BW_tmp(iX,iY);
        if part_ID>0
            BW_part=(BW_tmp==part_ID);
            part_pix=sum(sum(BW_part));
            part_int_val=sum(sum(HCO_stack_N(BW_part)))-part_pix*I_back_arr(iN);
            part_int_tmp=[part_int_tmp part_int_val];
        else
            part_int_tmp=[part_int_tmp 0];
        end
    end
    
    part_int=[part_int;part_int_tmp];
end

return

part_int_tmp=part_int;
for iP=1:size(part_int,1);part_int_tmp(iP,:)=part_int_tmp(iP,:)./median(part_int_tmp(iP,1:5));end
for iP=1:size(part_int,1);if and(max(part_int_tmp(iP,:))<20,min(part_int_tmp(iP,:))>0.05);plot(part_int_tmp(iP,:),'-');hold on;end;end

part_int_tmp=part_int;
part_stat=[];
for iP=1:size(part_int,1);part_int_tmp(iP,:)=part_int_tmp(iP,:)./median(part_int_tmp(iP,1:5));end
for iP=1:size(part_int,1);part_stat=[part_stat; median(part_int_tmp(iP,:)) std(part_int_tmp(iP,:))];end
figure
hist(part_stat(part_stat(:,1)<1.5,1),50)
figure
hist(part_stat(part_stat(:,2)<1,2),50)

part_int_tmp=part_int;
for iP=1:size(part_int,1);part_int_tmp(iP,:)=part_int_tmp(iP,:)./median(part_int_tmp(iP,1:5));end
for iP=1:size(part_int,1);if and(max(part_int_tmp(iP,:))<10,min(part_int_tmp(iP,:))>0.1);plot(polyval(pp,1./part_int_tmp(iP,:)*2.75),'-');hold on;end;end

for iP=1:size(part_int,1);if and(max(part_int_tmp2(iP,:))<10,min(part_int_tmp2(iP,:))>0.1);plot(part_int_tmp2(iP,:),'-');hold on;end;end
for iP=1:size(part_int,1);if and(max(part_int_tmp2(iP,:))<10,min(part_int_tmp2(iP,:))>0.1);plot(part_int_tmp2(iP,:),'-');pause;end;end
for iP=1:size(part_int,1);if and(max(part_int_tmp(iP,:))<10,min(part_int_tmp(iP,:))>0.1);plot(polyval(pp,1./part_int_tmp2(iP,:)*2.75),'-');pause;end;end
for iP=1:size(part_int,1);if and(max(part_int_tmp(iP,:))<10,min(part_int_tmp(iP,:))>0.1);plot(polyval(pp,1./part_int_tmp2(iP,:)*3.5),'-');pause;end;end

k_tmp=[];
for iP=1:size(part_int,1);if and(max(part_int_tmp(iP,:))<10,min(part_int_tmp(iP,:))>0.1);pp2=polyfit((5:50)*10,polyval(pp,1./part_int_tmp2(iP,5:50)*3.5),1);k_tmp=[k_tmp pp2(1)];;end;end
hist(k_tmp/3.5e-4,35)
for iP=1:size(part_int,1);if and(max(part_int_tmp(iP,:))<10,min(part_int_tmp(iP,:))>0.1);pp2=polyfit((5:40)*10,polyval(pp,1./part_int_tmp2(iP,5:40)*3.5),1);k_tmp=[k_tmp pp2(1)];;end;end
hist(k_tmp/3.75e-4,40)

k_tmp=[];
for iP=1:size(part_int,1);if and(max(part_int_tmp(iP,:))<10,min(part_int_tmp(iP,:))>0.1);pp2=polyfit((5:50)*10,polyval(pp,part_int_tmp2(iP,5:50)),1);k_tmp=[k_tmp pp2(1)];;end;end
hist(k_tmp(and(k_tmp>-1e-3,k_tmp<3.5e-3))/4e-4,29)
for iP=1:size(part_int,1);if and(max(part_int_tmp(iP,:))<10,min(part_int_tmp(iP,:))>0.1);plot(polyval(pp,part_int_tmp2(iP,:)),'-');pause;end;end


bleach_arr=[];
for iP=1:size(part_int,2);bleach_tmp=part_int_tmp(:,iP);bleach_arr=[bleach_arr median(bleach_tmp(isnan(bleach_tmp)==0))];end
for iP=1:size(part_int,2);bleach_tmp=part_int_tmp(:,iP);bleach_arr=[bleach_arr median(bleach_tmp(bleach_tmp>0))];end
bleach_arr(51:end)=bleach_arr(50);
bleach_arr(55:end)=bleach_arr(54);
part_int_tmp2=part_int;
for iP=1:size(part_int,1);part_int_tmp2(iP,:)=part_int_tmp2(iP,:)./median(part_int_tmp2(iP,1:5))./bleach_arr;end

part_int_back=part_int;
part_int=part_int_back(and(part_int_back(:,1)>700,part_int_back(:,1)<4000),:);

hist(power(part_int_back(:,1),1/3),50)

save(strcat(base_str,'.part_int.mat'),'part_int','-mat')

return

int_arr=[];
pH_arr=[];

for iP=1:size(ves_pos,1)
    tmptmp=tmp(round(ves_pos(iP,2)),round(ves_pos(iP,1)),:);tmptmp=reshape(tmptmp,size(tmptmp,3),1);    
    tmptmp=tmptmp-120;
    int_arr=[int_arr median(tmptmp(1:5))];
    if and(int_arr(end)>20,int_arr(end)<40)
        tmptmp=tmptmp./median(tmptmp(1:5));
        tmp_pH=polyval(pp,1./tmptmp*2.75);
%         if and(max(tmp_pH)<10,median(polyval(pp,1./tmptmp(50:60)*2.75))>6.5*0)
        plot(tmptmp,'-')
%           plot((0:102)*5,tmp_pH,'.')
%           tmp_pH_2=[];for iN=2:102;tmp_pH_2=[tmp_pH_2 mean(tmp_pH(iN-1:iN+1))];end
%           plot((1:101)*5,tmp_pH_2,'-')
%           pH_arr=[pH_arr tmptmp];
%           t_arr=[5.5+0.4 67.2+70.3 122+125.1 182.6+185.2 242.3+245.4 303.5+306 362.3+365.4 423.1+426.1 482+485.1 541.4+545 603.4+606.5 662.7+666.8 723.3+726.4 782.2+785.8]/2;
%           pH_dis_arr=[median(tmp_pH(1:11)) median(tmp_pH(12:18)) median(tmp_pH(19:25)) median(tmp_pH(26:31)) median(tmp_pH(32:38)) median(tmp_pH(39:44)) median(tmp_pH(45:51)) median(tmp_pH(52:58)) median(tmp_pH(59:65)) median(tmp_pH(66:73)) median(tmp_pH(74:80)) median(tmp_pH(81:89)) median(tmp_pH(90:96)) median(tmp_pH(97:104))];
%           plot(t_arr+iP/50,pH_dis_arr,'.')
          
          hold on
          pH_arr=[pH_arr median(polyval(pp,1./tmptmp(20:30)*2.75))];
%         end        
    end
end

% figure
% surf(log(tmp3'))
% shading interp
% view(0,90)

% figure
% tmptmp=tmp(109,564,:);tmptmp=reshape(tmptmp,size(tmptmp,3),1);
% plot(tmptmp-118,'.')
% hold on
% tmptmp=tmp(466,1691,:);tmptmp=reshape(tmptmp,size(tmptmp,3),1);
% plot(tmptmp-118,'.r')
clear;
iWind=1;

%%%
pp=[0.0206 -0.1963 0.9529 5.7134]; % 458
pp=[0.001 -0.027 0.247 -1.038 8.204];
pH_corr_2m0=0.94; % pH=6.45
r=75e-9;
pKs=7.2;

V=4/3*pi*power(r,3);
%%%

%%% 2 mM
base_str='210219-CcO-3zu1-EColipolar-cytc_0_5uM-004.nd2.1.nf-blur.tif';

tmp=load(strcat(base_str,'.part_int_mask.mat'));
part_int=tmp.part_int;

part_int_back=part_int;
part_int=part_int_back(and(part_int_back(:,1)>700/2,part_int_back(:,1)<2800*2),:);

% hist(power(part_int_back(:,1),1/3),50)

part_int_tmp=part_int;
for iP=1:size(part_int,1);part_int_tmp(iP,:)=part_int_tmp(iP,:)./median(part_int_tmp(iP,1:5));end
% for iP=1:size(part_int,1);if and(max(part_int_tmp(iP,:))<10,min(part_int_tmp(iP,:))>0.05);plot(part_int_tmp(iP,:),'-');hold on;end;end

% part_int_tmp2=part_int;
% for iP=1:size(part_int,1);part_int_tmp2(iP,:)=part_int_tmp2(iP,:)./median(part_int_tmp2(iP,1:5))./exp(polyval(pp_b,(1:size(part_int_tmp2,2))));end

% base_str='200805-CcO-3zu1-PCPE-pump-001.nd2.1.nf-blur.tif';
% 
% tmp=load(strcat(base_str,'.bleach.mat'));
% bleach_arr=tmp.bleach_arr(1:65);

bleach_arr=[];
for iP=1:size(part_int,2);bleach_tmp=part_int_tmp(:,iP);bleach_arr=[bleach_arr median(bleach_tmp(bleach_tmp>0))];end
bleach_arr_back=bleach_arr;
bleach_arr(56:end)=bleach_arr(55);
part_int_tmp2=part_int;
for iP=1:size(part_int,1);part_int_tmp2(iP,:)=part_int_tmp2(iP,:)./median(part_int_tmp2(iP,1:5))./bleach_arr*pH_corr_2m0;end

% figure
% plot(bleach_arr,'-')
% hold on
% plot(bleach_arr_back./bleach_arr,'-r')

k_tmp=[];
figure
for iP=1:size(part_int,1);
    if and(and(and(max(part_int_tmp(iP,:))<20,min(part_int_tmp(iP,:))>0.05),part_int_tmp(iP,50)>2),(part_int_tmp(iP,7)-1)./mean(part_int_tmp(iP,35:40)-1)<0.8);
        part_int_tmp3=polyval(pp,4./part_int_tmp2(iP,:));
%         plot((0:64)*10,(part_int_tmp2(iP,:)-1)./mean(part_int_tmp2(iP,40:45)-1),'-k');
        plot((0:64)*10,part_int_tmp3,'-k');
%         plot((0:64)*10,(part_int_tmp3(:)-polyval(pp,4))./mean(part_int_tmp3(35:40)-polyval(pp,4)),'-k');
%         hold on;
        close gcf
        pp2=polyfit((5:15)*10,polyval(pp,4./part_int_tmp2(iP,5:15)),1);k_tmp=[k_tmp pp2(1)];
    end;
end
axis([0 700 -0.5 2.5])
% axis([0 700 6.2 8.2])
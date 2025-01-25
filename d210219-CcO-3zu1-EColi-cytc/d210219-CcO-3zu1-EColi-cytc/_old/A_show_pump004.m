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
part_int=part_int_back(and(part_int_back(:,1)>700,part_int_back(:,1)<4000),:);
% part_int=part_int_back(and(part_int_back(:,1)>4000,part_int_back(:,1)<12000),:);

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
for iP=1:size(part_int,1);if and(max(part_int_tmp(iP,:))<20,min(part_int_tmp(iP,:))>0.1);pp2=polyfit((5:15)*10,polyval(pp,4./part_int_tmp2(iP,5:15)),1);k_tmp=[k_tmp pp2(1)];;end;end
figure
hsp=subplot(1,1,1);
[b,a]=hist(k_tmp,20);
bar(a,b./max(b)*100)
% plot(a,b./sum(b)*100,'-or','LineWidth',2)
% axis([-2 6 0 105])
set(hsp,'LineWidth',2)
set(hsp,'FontSize',24)
xlabel('slope [3.8 x 10^{-4}/s]')
ylabel('frequency [%]')

iWind=2;
std_arr=[];

figure
hsp=subplot(1,1,1);

for iP=1:size(part_int_tmp,1);
    if and(max(part_int_tmp(iP,:))<20,min(part_int_tmp(iP,:))>0.1);
        t_arr=[];
        I_arr=[];
        for iN=4+iWind:size(part_int_tmp(iP,:),2)-iWind
            t_arr=[t_arr mean((iN-iWind:iN+iWind))];
            I_arr=[I_arr mean(part_int_tmp2(iP,iN-iWind:iN+iWind))];
        end
%         plot(part_int_tmp(iP,:),'-');
        hold on;
        pp2=polyfit(t_arr(5:45)*10,polyval(pp,4./I_arr(5:45)),1);
        std_arr=[std_arr std(polyval(pp,4./I_arr(5:45))-polyval(pp2,t_arr(5:45)*10))];
        if std_arr(end)<0.1*100
            plot((t_arr(1:end-3)-iWind)*10,polyval(pp,4./I_arr(1:end-3)),'-','Color',0.6*[1 1 1],'LineWidth',1)
        end
%         pp2(1)
%         pause;
        hold on;
   end;
end
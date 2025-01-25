clear;
iWind=1;

%%% procedure
% 1. check parameters
% 2. check size selection
% 3. check bleaching
% 4. asign populations
% 5. check pH_out

%%% 2 mM
base_str='210219-CcO-3zu1-EColipolar-cytc_0_05uM-002.nd2.1.nf-blur.tif';
dt=1;

%%% vesicle parameters
r=90.3e-9;
lip_rat=49.3;
median_radius=10;
nig_int=(56:60);
pKs=7.1;
pHmax=7.5;
iWind=3;
% iWind=3; % for k(pH)
pp=polyfit([7.9 7.7 7.5 7.2 6.8 6.5 6.3 6.2],[6.8485 6.2565 5.0861 4.1940 2.0221 0.9955 0.8627 0.71],1);
%%% size selection
% - [186 201 207:221 223] most too fast, rest is inactive, 186 is slow
i_up=[186 201 207:221 223];
i_lower=350;
i_upper=700;
% - [185 195 197:220] most too fast, rest is inactive, 199 is slow
% fast [198 200 202 203 204 205 208 209 210 211 212 213 214 215 217 218 219 220] 18
% slow [199] 1
% chaotic [195 197 206 207] 4
% late [201 216] 2
i_up=[185 195 197:220]; % 0.1182
i_lower=700; % center
i_upper=1400;
% - [184:188 190:211] most too fast, rest is inactive, 184:185 is slow
% fast [186 187 188 190 191 192 193 194 195 196 197 198 199 200 201 202 203 205 206 207 208 209 210 211] 24
% slow [184 185] 2
% chaotic [204] 1
% late []
i_up=[184:188 190:211]; % 0.12809
i_lower=1400;
i_upper=2800;
% - [147 164 174 176 177 179:185 188:208] most too fast, rest is inactive,
% 164 is slow
% fast [] 31
% slow [] 0
% chaotic [177] 1
% late [174] 1
i_up=[147 164 176 179 180 181 182 183 184 185 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 208]; % fast
i_lower=2800;
i_upper=5600;

V=4/3*pi*power(r,3);
%%%

tmp=load(strcat('./processed/',base_str,'.part_int_mask.mat'));
part_int=tmp.part_int;
part_id=(1:size(part_int,1))';

part_int_back=part_int;
tmp_i_nig=median(part_int_back(:,nig_int),2);
tmp_i_base=median(part_int_back(:,1:5),2);
part_int=part_int_back(and(and(tmp_i_nig>i_lower,tmp_i_nig<i_upper),tmp_i_base>0),:);
part_id=part_id(and(and(tmp_i_nig>i_lower,tmp_i_nig<i_upper),tmp_i_base>0));

% hist(power(tmp,1/3),50)

part_int_tmp=part_int;
for iP=1:size(part_int,1);part_int_tmp(iP,:)=part_int_tmp(iP,:)./median(part_int_tmp(iP,nig_int));end
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
bleach_arr(nig_int(end)+1:end)=1;
% bleach_arr=1;
part_int_tmp2=part_int;
for iP=1:size(part_int,1);part_int_tmp2(iP,:)=part_int_tmp2(iP,:)./median(part_int_tmp2(iP,nig_int))./bleach_arr;end

% for iP=1:size(part_int,1);if and(max(part_int_tmp2(iP,:))<20,min(part_int_tmp2(iP,1:30))>0.05);plot(interp1(polyval(pp,(6:0.1:9.5)),(6:0.1:9.5),part_int_tmp2(iP,:)),'-');hold on;end;end

tmpX_all=[];tmpX_id_all=[];
for iP=1:size(part_int,1);if and(max(part_int_tmp2(iP,:))<20,min(part_int_tmp2(iP,1:30))>0.05);tmpX_all=[tmpX_all;interp1(polyval(pp,(6:0.1:9.5)),(6:0.1:9.5),part_int_tmp2(iP,:))];tmpX_id_all=[tmpX_id_all;part_id(iP)];end;end
tmpX_all=tmpX_all';
tmpX_id_all=tmpX_id_all';

% tmpX_all=tmpX_all(:,1:223);

sd=0*ones(size(tmpX_all,2),size(tmpX_all,2));
for iX=1:size(tmpX_all,2);
    for iY=iX+1:size(tmpX_all,2);
        tmptmp=tmpX_all(:,iX).*tmpX_all(:,iY);
        iarr=find(isnan(tmptmp)==0);
        sd(iX,iY)=mean(power(tmpX_all(iarr,iX)-tmpX_all(iarr,iY),2));
        sd(iY,iX)=sd(iX,iY);
    end;
end

tmptmp=median(sd,2);
[iy,ii]=sort(tmptmp);
% sdx=[];

sd=0*ones(size(tmpX_all,2),size(tmpX_all,2));
for iX=1:size(tmpX_all,2);
%     sdx=[sdx median(tmpX_all(iarr(1300:1400),ii(iX)))];
    for iY=1:size(tmpX_all,2);
        tmptmp=tmpX_all(:,ii(iX)).*tmpX_all(:,ii(iY));
        iarr=find(isnan(tmptmp)==0);
        sd(iX,iY)=mean(power(tmpX_all(iarr,ii(iX))-tmpX_all(iarr,ii(iY)),2));
    end;
end

% figure
% surf(log(sd))
% view(0,90)

sd2=0*ones(size(tmpX_all,2),size(tmpX_all,2));
for iX=1:size(tmpX_all,2);
    for iY=iX+1:size(tmpX_all,2);
        iarr=log(sd(:,iX))-log(sd(:,iY));
        iarr=iarr(isinf(iarr)==0);
        sd2(iX,iY)=mean(iarr.*iarr);
        sd2(iY,iX)=sd2(iX,iY);
    end;
end

% figure
% surf(sd2)
% view(0,90)

ii_tmp=ii;
ii2=[];
while length(ii)>0
    ii2=[ii(end) ii2];
    ii=ii(1:end-1);
    iarr1=find(ii_tmp==ii2(1));
    iarr2=find(sd2(iarr1,:)<0.2);
    for iX=1:length(iarr2)
       iarr3=find(ii==ii_tmp(iarr2(iX)));
       if length(iarr3)>0
           ii2=[ii(iarr3) ii2];
           ii(iarr3)=NaN;
           ii=ii(isnan(ii)==0);
       end
    end
end
ii=ii2';

sd=0*ones(size(tmpX_all,2),size(tmpX_all,2));
for iX=1:size(tmpX_all,2);
%     sdx=[sdx median(tmpX_all(iarr(1300:1400),ii(iX)))];
    for iY=1:size(tmpX_all,2);
        tmptmp=tmpX_all(:,ii(iX)).*tmpX_all(:,ii(iY));
        iarr=find(isnan(tmptmp)==0);
        sd(iX,iY)=mean(power(tmpX_all(iarr,ii(iX))-tmpX_all(iarr,ii(iY)),2));
    end;
end

sd2=0*ones(size(tmpX_all,2),size(tmpX_all,2));
for iX=1:size(tmpX_all,2);
    for iY=iX+1:size(tmpX_all,2);
        iarr=log(sd(:,iX))-log(sd(:,iY));
        iarr=iarr(isinf(iarr)==0);
        sd2(iX,iY)=mean(iarr.*iarr);
        sd2(iY,iX)=sd2(iX,iY);
    end;
end

% figure
% surf(sd2)
% view(0,90)

ii_tmp=ii;
ii2=[];
while length(ii)>0
    ii2=[ii(end) ii2];
    ii=ii(1:end-1);
    iarr1=find(ii_tmp==ii2(1));
    iarr2=find(sd2(iarr1,:)<0.5);
    for iX=1:length(iarr2)
       iarr3=find(ii==ii_tmp(iarr2(iX)));
       if length(iarr3)>0
           ii2=[ii(iarr3) ii2];
           ii(iarr3)=NaN;
           ii=ii(isnan(ii)==0);
       end
    end
end
ii=ii2';

sd=0*ones(size(tmpX_all,2),size(tmpX_all,2));
for iX=1:size(tmpX_all,2);
%     sdx=[sdx median(tmpX_all(iarr(1300:1400),ii(iX)))];
    for iY=1:size(tmpX_all,2);
        tmptmp=tmpX_all(:,ii(iX)).*tmpX_all(:,ii(iY));
        iarr=find(isnan(tmptmp)==0);
        sd(iX,iY)=mean(power(tmpX_all(iarr,ii(iX))-tmpX_all(iarr,ii(iY)),2));
    end;
end

sd2=0*ones(size(tmpX_all,2),size(tmpX_all,2));
for iX=1:size(tmpX_all,2);
    for iY=iX+1:size(tmpX_all,2);
        iarr=log(sd(:,iX))-log(sd(:,iY));
        iarr=iarr(isinf(iarr)==0);
        sd2(iX,iY)=mean(iarr.*iarr);
        sd2(iY,iX)=sd2(iX,iY);
    end;
end

% figure
% surf(sd2)
% view(0,90)

% figure;plot(bleach_arr_back)

% figure;for iX=[147 164 174 176 177 179:185 188:208];plot((1:size(tmpX_all,1))*dt,tmpX_all(:,ii(iX)),'-','Tag',num2str(tmpX_id_all(ii(iX))));hold on;end;axis([0 70*dt 6 9.5])

%%%

% figure
k_H_tmp=[];

for iP=1:length(i_up)
    r_act=power(median(part_int_back(tmpX_id_all(ii(i_up(iP))),nig_int)),1/3)*r./median_radius;
    N_lum=4/3*pi*power(r_act,3)*2.4e-3*6e26;
    N_lip=4*pi*power(r_act,2)./0.64e-18/lip_rat;
    buffer_conc=(N_lum+N_lip)./N_lum.*2.4e-3; % M
    V=4/3*pi*power(r_act,3);
    
    log_chi=(-1:0.01:1)*2;
    chi=power(10,log_chi);
    prot_buffer=chi./(1+chi);
    pH=pKs-log_chi;
    free_H=power(10,-pH);               % conc. unbound H+h
    buffer_H=buffer_conc*prot_buffer;   % conc. protonated buffer
    total_H=buffer_H+free_H;            % conc. added H+ to achieve desired pH value
    
%     figure
%     hsp=subplot(1,1,1);

    t_arr=[];
    I_arr=[];
    tmpX=(1:size(tmpX_all,1))*dt;
    tmpY=tmpX_all(:,ii(i_up(iP)))';
    for iN=4+iWind:max(nig_int)-iWind
        t_arr=[t_arr mean(tmpX(iN-iWind:iN+iWind))];
        I_arr=[I_arr mean(tmpY(iN-iWind:iN+iWind))];
    end
    iarr=find(I_arr(1:end-10)<pHmax);
    pp2=polyfit(t_arr(iarr),interp1(-log10(free_H),(buffer_H+free_H)*6e26*V,I_arr(iarr)),1);k_H_tmp=[k_H_tmp pp2(1)];

%     [i_up(iP) pp2(1)]
%     subplot(2,1,1)
%     plot(t_arr,I_arr,'-','Tag',num2str(tmpX_id_all(ii(iX))));
%     subplot(2,1,2)

%     plot(t_arr-iWind,interp1(-log10(free_H),(buffer_H+free_H)*6e26*V,I_arr),'-')

%     plot(t_arr,interp1(-log10(free_H),(buffer_H+free_H)*6e26*V,I_arr),'-','Color',0.6*[1 1 1])
%     hold on
%     plot(t_arr(iarr),polyval(pp2,t_arr(iarr)),'-r')
    
    tmp_diff=-diff(interp1(-log10(free_H),(buffer_H+free_H)*6e26*V,I_arr(1:end-10)))./dt;
    semilogy(I_arr(1:end-11),tmp_diff,'o','Color',0.6*[1 1 1],'LineWidth',1.5)
%     semilogy(I_arr(1:end-11),tmp_diff,'ob','LineWidth',1.5)
%     semilogy(I_arr(1:end-11),tmp_diff,'-')
%     semilogy(I_arr(3:end-10)-6.4,tmp_diff,'o')
    hold on
    
    close gcf
end
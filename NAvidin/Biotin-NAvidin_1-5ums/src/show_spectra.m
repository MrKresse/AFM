
function show_spectra(path_processed, noise_str)
NoPeaks=1;
PeakArr=[0 4e4];
L=349.69e-6;   % 350µm lang, copy C096
w=32.43e-6;    %  35µm breit, copy C096
invols=110e-9;
%

%NData=64000;
kB=1.3806505e-23;
T=300;

%E=150e9;    % Elastizitätsmodul, Si
E=169e9;    % Elastizitätsmodul, Si
rho=2.33e3; % Dichte Beam, Si
%rho_fl=1.18; % Dichte Umgebung, Luft
%eta=1.86e-5;   % Viskosität Umgebung, Luft
rho_fl = 997; % Dichte Wasser
eta =1e-3; % Viskosität wasser


file_str=strcat(path_processed,noise_str);
fid = fopen(file_str,'r','s');
noise_spec=fread(fid,500000,'single');
fclose(fid);

%noise_spec=noise_spec-1e-10;            %% AD-noise, Taken from VI
noise_spec=noise_spec*invols*invols;

if (1==0)
    semilogy((1:500000),noise_spec,'.');
    return;
end

for nJ=1:NoPeaks
  %Mitte bestimmen
  [k,l]=max(noise_spec(PeakArr(nJ,1):PeakArr(nJ,2)));
  l=l+PeakArr(nJ,1)-1;
  pp=polyfit((l-300:l+300)'-1,noise_spec(l-300:l+300),2);
  x0=-pp(2)./(2*pp(1));
  s=1./(pp(3)-pp(2)*pp(2)/4/pp(1));
  
  if nJ==1
    %Mitte bestimmen
    y0=polyval(pp,x0);
    PeakArrTemp=noise_spec(PeakArr(nJ,1):PeakArr(nJ,2))./y0;

    std_arr=[];
    Q_arr=[];

    for Q=(1:.02:40)
      x_a=1./(power(4*pi*pi*x0.*x0-4*pi*pi*((PeakArr(nJ,1):PeakArr(nJ,2))-1).*((PeakArr(nJ,1):PeakArr(nJ,2))-1),2)+power(4*pi*pi*x0.*((PeakArr(nJ,1):PeakArr(nJ,2))-1)/Q,2));
      x_a=x_a./max(x_a);
      std_arr=[std_arr sum(abs(PeakArrTemp'-x_a))];
      Q_arr=[Q_arr Q];
    end

    [y,m]=min(std_arr);
    Q=Q_arr(m);
    x_a=1./(power(4*pi*pi*x0.*x0-4*pi*pi*((PeakArr(nJ,1):PeakArr(nJ,2))-1).*((PeakArr(nJ,1):PeakArr(nJ,2))-1),2)+power(4*pi*pi*x0.*((PeakArr(nJ,1):PeakArr(nJ,2))-1)/Q,2));
    x_a=x_a./max(x_a);
  
    Re=rho_fl.*(2*pi*x0).*w.*w./4/eta;
    tau=log10(Re);
    ReOmega=polyval([.00069085 -0.0035117 0.044055 -0.12886 0.46842 -0.48274 0.91324],tau)./polyval([.00069085 -0.0035862 0.045155 -0.13444 0.4869 -0.56964 1],tau);
    ImOmega=polyval([-0.000044510 + 0.000064577 -0.00010961 0.016294 -0.029256 -0.024134],tau)./polyval([0.00286361 -0.014369 0.079156 -0.18357 0.55182 -0.59702 1],tau);
    Omega=ReOmega+i*ImOmega;

    Gamma_Circ=1+4*i*besselk(1,-i*sqrt(i*Re))./(sqrt(i*Re).*besselk(0,-i*sqrt(i*Re)));
    Gamma=Omega.*Gamma_Circ;

    figure
    hsp=subplot(1,1,1);
    semilogy(((PeakArr(nJ,1):PeakArr(nJ,2))-1)/1e3,PeakArrTemp*y0,'.','MarkerSize',10)
    hold on
    semilogy(((PeakArr(nJ,1):PeakArr(nJ,2))-1)/1e3,x_a*y0,'r','LineWidth',2)
    set(hsp,'FontSize',16);
    set(hsp,'LineWidth',2)
    xlabel('frequency [kHz]');
    ylabel('noise density [m^2/Hz]');
    hlg=legend('measurement',strcat('lorentzian fit, k_{Sader}=',num2str(0.1906*L*Q*power(2*pi*x0*w,2)*rho_fl*imag(Gamma)),', k_{Cleveland}=',num2str(2*w*power(pi*L*x0,3)*sqrt(power(rho,3)./E)),', Q=',num2str(Q)));
    set(hlg,'FontSize',13);
    set(hlg,'Box','off');
    set(gcf,'PaperOrientation','portrait');
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'PaperPosition',[0.634518 0.621759 19.715 13.9237]);    
    set(gcf,'PaperSize',[20.984 15.18]);
    
    filestr=strcat('processed/',noise_str,'.',num2str(nJ),'.Q.png');
    print('-dpng',filestr);
    filestr=strcat('processed/',noise_str,'.',num2str(nJ),'.Q.pdf');
    print('-dpdf',filestr);
    hold off
  end  
  
  marr=[];
  iarr=(1:50);
  for j=1:50;
    marr=[marr mean(abs(noise_spec(l-400:l+400)'-s./(s.*s+power(10,j)*power((l-400:l+400)-1-x0,2))))];
  end;

  [k,i]=min(abs(marr));
  base1=iarr(i);

  marr=[];
  iarr=(0:0.005:10);
  for j=0:0.005:10;
    marr=[marr mean(abs(noise_spec(l-400:l+400)'-s./(s.*s+j*power(10,base1)*power((l-400:l+400)-1-x0,2))))];
  end;
  [k,i]=min(abs(marr));
  base2=iarr(i);

  a=base2*power(10,base1);

  figure
  hsp=subplot(1,1,1);
  semilogy(((PeakArr(nJ,1):PeakArr(nJ,2))-1)/1e3,noise_spec(PeakArr(nJ,1):PeakArr(nJ,2)),'.','MarkerSize',10)
  hold on
  semilogy(((PeakArr(nJ,1):PeakArr(nJ,2))-1)/1e3,s./(s.*s+a*power(((PeakArr(nJ,1):PeakArr(nJ,2))-1)-x0,2)),'r','LineWidth',2)
  set(hsp,'FontSize',16);
  set(hsp,'LineWidth',2)
  xlabel('frequency [kHz]');
  ylabel('noise density [m^2/Hz]');
  if nJ==1
    hlg=legend('measurement',strcat('lorentzian fit, k=',num2str(0.8174*kB*T/(pi/sqrt(a)))),'jamming limit');
  elseif nJ==2
    hlg=legend('measurement',strcat('lorentzian fit, k=',num2str(0.2511*kB*T/(pi/sqrt(a)))),'jamming limit');
  elseif nJ==3
    hlg=legend('measurement',strcat('lorentzian fit, k=',num2str(0.0863*kB*T/(pi/sqrt(a)))),'jamming limit');
  end      
  set(hlg,'FontSize',13);
  set(hlg,'Box','off');
  set(gcf,'PaperOrientation','portrait');
  set(gcf,'PaperPositionMode','manual');
  set(gcf,'PaperPosition',[0.634518 0.621759 19.715 13.9237]);    
  set(gcf,'PaperSize',[20.984 15.18]);
  
  filestr=strcat('processed/',noise_str,'.',num2str(nJ),'.L.png');
  print('-dpng',filestr);
  filestr=strcat('processed/',noise_str,'.',num2str(nJ),'.L.pdf');
  print('-dpdf',filestr);
  hold off
end

%close(gcf);

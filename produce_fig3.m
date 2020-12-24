clear all;

N=100;
kn=linspace(0,5,N);
for j=1:N
  disp(j);
  
  %% case1-plots in blue. 
  g1d=1; g1e=1; g2d=1; g2e=1;
  [Ef1, Ef2, vg1, vg2, nc]=thermalshg(kn(j),g1d,g1e,g2d,g2e);
  Efb1(j)=Ef1; Efb2(j)=Ef2; vgb1(j)=vg1; vgb2(j)=vg2; ncb(j)=nc;

  %% case2-yellow plots
  g1d=1; g1e=0.1; g2d=1; g2e=1;
  [Ef1, Ef2, vg1, vg2, nc]=thermalshg(kn(j),g1d,g1e,g2d,g2e);
  Efy1(j)=Ef1; Efy2(j)=Ef2; vgy1(j)=vg1; vgy2(j)=vg2; ncy(j)=nc;

  %% case3-green plots
  g1d=1; g1e=1; g2d=0.1; g2e=1;
  [Ef1, Ef2, vg1, vg2, nc]=thermalshg(kn(j),g1d,g1e,g2d,g2e);
  Efg1(j)=Ef1; Efg2(j)=Ef2; vgg1(j)=vg1; vgg2(j)=vg2; ncg(j)=nc;

  %% case4-red plots
  g1d=1; g1e=0.1; g2d=0.1; g2e=1;
  [Ef1, Ef2, vg1, vg2, nc]=thermalshg(kn(j),g1d,g1e,g2d,g2e);
  Efr1(j)=Ef1; Efr2(j)=Ef2; vgr1(j)=vg1; vgr2(j)=vg2; ncr(j)=nc;
  
end

figure(1);  % enhancement factor at frequency w2
plot(kn,Efb2,'b',kn,Efy2,'y',kn,Efg2,'g',kn,Efr2,'r','LineWidth',2);
xlabel('\kappa/\gamma'); ylabel('F_2'); title('Emission enhancement at \omega_2');

figure(2);
plot(kn,Efb1,'b',kn,Efy1,'y',kn,Efg1,'g',kn,Efr1,'r','LineWidth',2);
xlabel('\kappa/\gamma'); ylabel('F_1'); title('Emission enhancement at \omega_1');

figure(3);  % g2 value at w2 
plot(kn,vgb2,'b',kn,vgy2,'y',kn,vgg2,'g',kn,vgr2,'r','LineWidth',2);
xlabel('\kappa/\gamma'); ylabel('g_2^{(0)} at \omega_2'); title('Statistics of emission at \omega_2');

figure(4);  % g2 value at w2 
plot(kn,vgb1,'b',kn,vgy1,'y',kn,vgg1,'g',kn,vgr1,'r','LineWidth',2);
xlabel('\kappa/\gamma'); ylabel('g_2^{(0)} at \omega_1'); title('Statistics of emission at \omega_1');

figure(5);  % intensity correlations
plot(kn,ncb,'b',kn,ncy,'y',kn,ncg,'g',kn,ncr,'r','LineWidth',2);
xlabel('\kappa/\gamma'); ylabel('C'); title('Biphoton intensity correlations');





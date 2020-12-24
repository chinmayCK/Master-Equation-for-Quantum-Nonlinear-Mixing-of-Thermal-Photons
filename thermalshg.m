function[Ef1, Ef2, vg1, vg2, nc]=thermalshg(kn,g1d,g1e,g2d,g2e)

  c=2.997928e8; 
  w1=2*pi*c/(12e-6);  % SPP frequency 
  w2=2*w1;  % resonator frequency 
  Td=600;  % device temperaure

  n1=navg(w1,Td);  % avg. photon number at w1 for temperature =Td
  n2=navg(w2,Td);
  nt=n1+n2;
      
  % Dissipative intrinsic bath 
  Ga1=g1d; % resonator dissipative decay
  Ga2=g2d; % resonator dissipative decay

  % radiative environment
  Te=0.01;  % external environment temperature
  Ga1r=g1e; % radiative decay 
  Ga2r=g2e; % radiative decay 
  n1r=navg(w1,Te);
  n2r=navg(w2,Te);

  % in the linear regime, the effective photon number is decay-rate-weighted avg. 
  n2avg=(n2*Ga2+n2r*Ga2r)/(Ga2+Ga2r); 
  n1avg=(n1*Ga1+n1r*Ga1r)/(Ga1+Ga1r);
  
  N=10;

  for j=1:N
    % annihilation operator 
    va(j)=sqrt(j);
  end
  for j=1:(N+1)
    n=j-1;
    pn1(j)=1/(n1+1)*(n1/(1+n1))^n;  % this corresponds to pa linear
    pn2(j)=1/(n2+1)*(n2/(1+n2))^n;  % this corresponds to pa linear 
  end
  Da1m=sparse(diag(va,1));   Da1p=sparse(diag(va,-1));
  Da2m=Da1m;   Da2p=Da1p;
  Ia1=speye(N+1); Ia2=Ia1;

  %pa1=diag(pn1);
  %pa2=diag(pn2);

  %v_ar=kron(pa1(:),pa2(:));
  I_ar=kron(Ia1(:),Ia2(:));

  % Terms corresponding to coupling between the modes in master equation
  T1=1i*kn*kron(kron(Ia1,(Da1p*Da1p).'),kron(Ia2,Da2m.'));
  T2=1i*kn*kron(kron(Ia1,(Da1m*Da1m).'),kron(Ia2,Da2p.'));
  T3=-1i*kn*kron(kron(Da1p*Da1p,Ia1),kron(Da2m,Ia2));
  T4=-1i*kn*kron(kron(Da1m*Da1m,Ia1),kron(Da2p,Ia2));

  % Terms coresponding to linear Lindbladian [intrinsic disispation of a1]
  Lda1=Ga1*(n1+1)*(2*kron(kron(Da1m,Da1p.'),kron(Ia2,Ia2))...
		  -kron(kron(Da1p*Da1m,Ia1),kron(Ia2,Ia2))...
		  -kron(kron(Ia1,(Da1p*Da1m).'),kron(Ia2,Ia2)))...
      +Ga1*n1*(2*kron(kron(Da1p,Da1m.'),kron(Ia2,Ia2))...
		  -kron(kron(Da1m*Da1p,Ia1),kron(Ia2,Ia2))...
		  -kron(kron(Ia1,(Da1m*Da1p).'),kron(Ia2,Ia2)));

  % intrinsic dissipation of a2 
  Lda2=Ga2*(n2+1)*(2*kron(kron(Ia1,Ia1),kron(Da2m,Da2p.'))...
		  -kron(kron(Ia1,Ia1),kron(Da2p*Da2m,Ia2))...
		  -kron(kron(Ia1,Ia1),kron(Ia2,(Da2p*Da2m).')))...
      +Ga2*n2*(2*kron(kron(Ia1,Ia1),kron(Da2p,Da2m.'))...
		  -kron(kron(Ia1,Ia1),kron(Da2m*Da2p,Ia2))...
		  -kron(kron(Ia1,Ia1),kron(Ia2,(Da2m*Da2p).')));

  % radiative decay of a1
  Lra1=Ga1r*(n1r+1)*(2*kron(kron(Da1m,Da1p.'),kron(Ia2,Ia2))...
		  -kron(kron(Da1p*Da1m,Ia1),kron(Ia2,Ia2))...
		  -kron(kron(Ia1,(Da1p*Da1m).'),kron(Ia2,Ia2)))...
       +Ga1r*n1r*(2*kron(kron(Da1p,Da1m.'),kron(Ia2,Ia2))...
		  -kron(kron(Da1m*Da1p,Ia1),kron(Ia2,Ia2))...
		  -kron(kron(Ia1,(Da1m*Da1p).'),kron(Ia2,Ia2)));

  % radiative decay of a2
  Lra2=Ga2r*(n2r+1)*(2*kron(kron(Ia1,Ia1),kron(Da2m,Da2p.'))...
		  -kron(kron(Ia1,Ia1),kron(Da2p*Da2m,Ia2))...
		  -kron(kron(Ia1,Ia1),kron(Ia2,(Da2p*Da2m).')))...
      +Ga2r*n2r*(2*kron(kron(Ia1,Ia1),kron(Da2p,Da2m.'))...
		  -kron(kron(Ia1,Ia1),kron(Da2m*Da2p,Ia2))...
		  -kron(kron(Ia1,Ia1),kron(Ia2,(Da2m*Da2p).')));
  
  L=1*(T1+T2+T3+T4)+Lda1+Lda2+Lra1+Lra2; 
  [m,n]=size(L);
  L(m+1,:)=I_ar.';
  B=zeros(m+1,1); B(end)=1;  

  S=L\B;
  
  I_ar.'*S;   % check if the trace is really one 

  Dna1=Da1p*Da1m;   % for <a1p*a1m> calculation 
  Dna2=Da2p*Da2m;   % for <a2p*a2m> calculation 
  Ds1=Da1p*Da1p*Da1m*Da1m;   % for g2 at w1 
  Ds2=Da2p*Da2p*Da2m*Da2m;   % for g2 at w2 
  
  va1=kron(Dna1(:),Ia2(:));  % for <a1p*a1m> calculation 
  va2=kron(Ia1(:),Dna2(:));  % for <a2p*a2m> calculation 

  vc=kron(Dna1(:),Dna2(:));  % intensity corrleations 
  vg1=kron(Ds1(:),Ia2(:));  % for g2 calculation
  vg2=kron(Ia1(:),Ds2(:));  % for g2 calculation 
  
  nn1=va1.'*S;  % <a1p*a1m>
  nn2=va2.'*S;   % <a2p*a2m>
  nc=real(vc.'*S/(nn1*nn2)-1);     % intensity correlations strength 

  vg1=real(vg1.'*S/(nn1^2));  % g2 at w1
  vg2=real(vg2.'*S/(nn2^2));  % g2 at w2
  
  % interaction energy calculation
  D1p=Da1p*Da1p;  D1m=Da1m*Da1m;
  v1=kron(D1p(:),Da2m(:));
  v2=kron(D1m(:),Da2p(:));
  vi1=v1.'*S;  % <a1p*a1p*a2m> interaction 
  vi2=v2.'*S;  % <a1m*a1m*a2p> interaction  
  
  %rho=new_reshape(S);
  %spy(rho)
  %entl=Negativity(rho)

  Ef2=real(nn2/n2avg);   % enhancement factor of emission at w2
  Ef1=real(nn1/n1avg);   % enhancement/suppresion factor of emission at w1
% the quantities nn1,nn2 are treated as complex even though imag part is zero. Due to numerical methods in solving
% master equation.
  
  Ef=nn1/n1r;

  return;

function[A]=navg(w,T)
  hbar=1.05457e-34; 
  kb=1.38e-23; 
  A=1./(exp(abs(hbar*w/(kb*T)))-1);
return;

function[A]=abs2(a)
  A=a.*conj(a);
return;

function[rho]=new_reshape(vrho,dim1,dim2)
   if nargin <2
     dim1 = sqrt(sqrt(length(vrho)));
     dim2=dim1;
   end
   rho = sparse(dim1*dim2,dim1*dim2);
   for n=1:dim1^2
     idx1=dim2^2*(n-1)+1;
     idx2=dim2^2*n;
     vec = vrho(idx1:idx2);
     rhon = reshape(vec,[dim2,dim2]);
     
     vn=zeros(dim1^2,1); vn(n)=1;
     idm1=reshape(vn,[dim1,dim1]);
     rho=rho+kron(sparse(idm1),rhon);
   end 
  return

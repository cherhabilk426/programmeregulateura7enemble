
clc;clear all;close all;

%%%%%%%%%%%%%%%%%%%%% Paramètres de la machine 

rs=4.85;rr=3.805;ls=0.274;lr=0.274;lm=0.258;
j=0.031;p=2;f=0.008;
tr=lr/rr;
ts=ls/rs;
s=1-lm*lm/(lr*ls);

%*********************************************************************************************************

f1=0;f2=0;f3=0;f4=0;f5=0;

h=0.0001;

%%%%%%%%%% Paramètres de l'éstimateur

k_est=2;
Te_est=k_est*h;
phdrestim=0.1;
phqrestim=0.1; 
phrestim=0.001;

n1est=0;n2est=0;

%************* Paramètres des régulateurs flous ******************************************************
n1=0;n2=0;
nn1=0;nn2=0;

iqsmax=13;
idsmax=5;

k_regw=10;
k_regphi=5;

ew0=0;
ephr0=0;

idsref=0;
iqsref=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tf=1.5;

%----------------------------------
j_w=0;
j_phr=0;

k=1;
t=k*h;
x=[0;0;0;0;0;0];

while t<=tf;
    
   %cr=0;
   if (t > 0.5 & t < 1) cr=10; else cr=0; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wsg=lm*x(4)/(tr*phrestim);
ws=x(5)+wsg;
  
ids=x(3);
iqs=x(4);

%%%%%%%%% Estimateur de flux *********************************************************************************

if n1est==n2est
   
phdrestim=phdrestim+Te_est*(lm*ids+tr*wsg*phqrestim-phdrestim)/tr;
phqrestim=phqrestim+Te_est*(lm*iqs-tr*wsg*phdrestim-phqrestim)/tr;

phrestim=sqrt(phdrestim*phdrestim+phqrestim*phqrestim);
phrestim0=phrestim;

end;

%************** Régulation floue de la vitesse ********************************************************

  Nref=1000;
  %if t < 0.7, Nref=1000; else Nref=-1000; end;

N=30*x(5)/(p*pi);

if n1==n2
    
%------------------------
j_w=j_w+1;
%-----------------------
    
ew=(Nref-N)/abs(Nref);
vew=(ew-ew0)*50;

ew0=ew;

  [iqsref0,diqsref]=Floue7_w(iqsref,ew,vew);

  diqsref0(j_w)=diqsref;
  ew00(j_w)=ew;
  vew00(j_w)=vew;

end;

if abs(iqsref0)> iqsmax,
   
   iqsref=iqsmax*sign(iqsref0);
   
else
   
  iqsref=iqsref0;

end;

%************** Régulation floue du flux *****************************************************************

   phrref=1;

if nn1==nn2
    
%--------------
j_phr=j_phr+1;
%-------------
   
ephr=(phrref-phrestim)/1;
vephr=(ephr-ephr0)*40;
ephr0=ephr;

[idsref0,didsref]=Floue7_phr(idsref,ephr,vephr);

didsref0(j_phr)=didsref;
ephr00(j_phr)=ephr;
vephr00(j_phr)=vephr;

end;

if abs(idsref0)> idsmax,
   
   idsref=idsmax*sign(idsref0);
   
else
   
   idsref=idsref0;
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thetas=x(6);

iasref=sqrt(2/5)*(idsref*cos(thetas)-iqsref*sin(thetas));
ibsref=sqrt(2/5)*(idsref*cos(thetas+2*pi/5)-iqsref*sin(thetas+2*pi/5));
icsref=sqrt(2/5)*(idsref*cos(thetas+4*pi/5)-iqsref*sin(thetas+4*pi/5));
ifsref=sqrt(2/5)*(idsref*cos(thetas+6*pi/5)-iqsref*sin(thetas+6*pi/5));
iesref=sqrt(2/5)*(idsref*cos(thetas+8*pi/5)-iqsref*sin(thetas+8*pi/5));


ias=sqrt(2/5)*(x(3)*cos(thetas)-x(4)*sin(thetas));
ibs=sqrt(2/5)*(x(3)*cos(thetas+2*pi/5)-x(4)*sin(thetas+2*pi/5));
ics=sqrt(2/5)*(x(3)*cos(thetas+4*pi/5)-x(4)*sin(thetas+4*pi/5));
ifs=sqrt(2/5)*(x(3)*cos(thetas+6*pi/5)-x(4)*sin(thetas+6*pi/5));
ies=sqrt(2/5)*(x(3)*cos(thetas+8*pi/5)-x(4)*sin(thetas+8*pi/5));

tensionref=ondhyst(ias,ibs,ics,ifs,ies,   iasref,ibsref,icsref,ifsref,iesref,   f1,f2,f3,f4,f5);

vasref=tensionref(1);
vbsref=tensionref(2);
vcsref=tensionref(3);
vfsref=tensionref(4);
vesref=tensionref(5);

f1=tensionref(6);  %??  4
f2=tensionref(7);
f3=tensionref(8);
f4=tensionref(9);
f5=tensionref(10);

vds= sqrt(2/5)*(vasref*cos(thetas)+vbsref*cos(thetas+2*pi/5)+vcsref*cos(thetas+4*pi/5) +vfsref*cos(thetas+6*pi/5)+vesref*cos(thetas+8*pi/5));
vqs=-sqrt(2/5)*(vasref*sin(thetas)+vbsref*sin(thetas+2*pi/5)+vcsref*sin(thetas+4*pi/5) +vfsref*sin(thetas+6*pi/5)+vesref*sin(thetas+8*pi/5));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1=x;                         k1=derives(x1,vds,vqs,ws,cr);
x2=x1+0.5*h*k1;               k2=derives(x2,vds,vqs,ws,cr);
x3=x1+0.5*h*k2;               k3=derives(x3,vds,vqs,ws,cr); 
x4=x1+h*k3;                   k4=derives(x4,vds,vqs,ws,cr);

x=x1+h*(k1+2*k2+2*k3+k4)/6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phds(k)=x(1);
phqs(k)=x(2);
ids(k)=x(3);
iqs(k)=x(4);
wm(k)=x(5);
thetas=x(6);

N0(k)=N;
Nref0(k)=Nref;

cem(k)=p*(iqs(k)*phds(k)-ids(k)*phqs(k));
phdr(k)=lr*phds(k)/lm-s*ls*lr*ids(k)/lm;
phqr(k)=lr*phqs(k)/lm-s*ls*lr*iqs(k)/lm;
phr(k)=sqrt(phdr(k)*phdr(k)+phqr(k)*phqr(k));
iass(k)=sqrt(2/3)*(ids(k)*cos(thetas)-iqs(k)*sin(thetas));

phdrestim0(k)=phdrestim;
phqrestim0(k)=phqrestim;

%******************************************************************************************************************

iqsref00(k)=iqsref;

idsref00(k)=idsref;

n1=k/k_regw;
n2=fix(n1);

nn1=k/k_regphi;
nn2=fix(nn1);

n1est=k/k_est;
n2est=fix(n1est);

tk(k)=t; home, t
k=k+1; t=t+h;

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xx1=[ephr00;vephr00]';
xx2=[ew00;vew00]';

didsref00=didsref0';
diqsref00=diqsref0';

save datad xx1 didsref00;
save dataq xx2 diqsref00;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot(3,2,1),plot(tk,N0,tk,Nref0),grid, ylabel('n(tr/m)');
axis([0 tf 0 1100])
%axis([0 tf -1100 1100])
subplot(3,2,2),plot(tk,phdr,tk,phdrestim0),grid, ylabel('phdr(wb)'); 
axis([0 tf 0 1.1])
subplot(3,2,3),plot(tk,cem),grid, ylabel('cem(Nm)'); 
subplot(3,2,4),plot(tk,phqr,tk,phqrestim0),grid, ylabel('phqr(wb)'); 
subplot(3,2,5),plot(tk,iass),grid,xlabel('t(s)'); ylabel('ias(A)'); 
subplot(3,2,6),plot(tk,iqsref00,tk,idsref00),grid,xlabel('t(s)'); ylabel('idsref et iqsref(A)');

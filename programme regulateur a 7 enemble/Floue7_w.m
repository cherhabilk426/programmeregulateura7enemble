function [iqsref,diqsref]=Floue7_w(iqsref,ew,vew);

a1=0.25;a2=0.5;a3=1;

% fonction mu (NG,NM,NP,ZE,PP,PM,PG) pour l'erreur

muewNG=max(min(1,-(ew+a2)/(a3-a2)),0);
muewNM=max(min((ew+a3)/(a3-a2),-(ew+a1)/(a2-a1)),0);
muewNP=max(min((ew+a2)/(a2-a1),-ew/a1),0);
muewZE=max(min(1+ew/a1,1-ew/a1),0);
muewPP=max(min(ew/a1,(a2-ew)/(a2-a1)),0);
muewPM=max(min((ew-a1)/(a2-a1),(a3-ew)/(a3-a2)),0);
muewPG=max(min((ew-a2)/(a3-a2),1),0);

%fonction mu (NG,NM,NP,ZE,PP,PM,PG) pour la dérivée

muvewNG=max(min(1,-(vew+a2)/(a3-a2)),0);
muvewNM=max(min((vew+a3)/(a3-a2),-(vew+a1)/(a2-a1)),0);
muvewNP=max(min((vew+a2)/(a2-a1),-vew/a1),0);
muvewZE=max(min(1+vew/a1,1-vew/a1),0);
muvewPP=max(min(vew/a1,(a2-vew)/(a2-a1)),0);
muvewPM=max(min((vew-a1)/(a2-a1),(a3-vew)/(a3-a2)),0);
muvewPG=max(min((vew-a2)/(a3-a2),1),0);

% defuzzification : max-min
muqNG=max([min(muewNG,muvewNG);min(muewNG,muvewNM);min(muewNG,muvewNP);
           min(muewNG,muvewZE);min(muewNM,muvewNG);min(muewNM,muvewNM);
           min(muewNM,muvewNP);min(muewNP,muvewNG);min(muewNP,muvewNM);min(muewZE,muvewNG)]);
muqNM=max([min(muewNG,muvewPP);min(muewNM,muvewZE);min(muewNP,muvewNP);
           min(muewZE,muvewNM);min(muewPP,muvewNG)]);
muqNP=max([min(muewNG,muvewPM);min(muewNM,muvewPP);min(muewNP,muvewZE);
           min(muewZE,muvewNP);min(muewPP,muvewNM);min(muewPM,muvewNG)]);
muqZE=max([min(muewNG,muvewPG);min(muewNM,muvewPM);min(muewNP,muvewPP);
           min(muewZE,muvewZE);min(muewPP,muvewNP);min(muewPM,muvewNM);min(muewPG,muvewNG)]);
muqPP=max([min(muewNM,muvewPG);min(muewNP,muvewPM);min(muewZE,muvewPP);
           min(muewPP,muvewZE);min(muewPM,muvewNP);min(muewPG,muvewNM)]);
muqPM=max([min(muewNP,muvewPG);min(muewZE,muvewPM);min(muewPP,muvewPP);
           min(muewPM,muvewZE);min(muewPG,muvewNP)]);
muqPG=max([min(muewZE,muvewPG);min(muewPP,muvewPG);min(muewPP,muvewPM);
           min(muewPM,muvewPG);min(muewPM,muvewPM);min(muewPM,muvewPP);
           min(muewPG,muvewPG);min(muewPG,muvewPM);min(muewPG,muvewPP);min(muewPG,muvewZE)]);
       

som_muq=muqNG+muqNM+muqNP+muqZE+muqPP+muqPM+muqPG;

%Iq=13;Guq=0.95;

Iq=1;

s_PG=1*Iq;
s_PM=0.5*Iq;
s_PP=0.25*Iq;
s_ZE=0*Iq;
s_NP=-0.25*Iq;
s_NM=-0.5*Iq;
s_NG=-1*Iq;

diqsref=(s_PG*muqPG+s_PM*muqPM+s_PP*muqPP+s_ZE*muqZE+s_NP*muqNP+s_NM*muqNM+s_NG*muqNG)/som_muq;

Guq=6.5;
iqsref=iqsref+Guq*diqsref;

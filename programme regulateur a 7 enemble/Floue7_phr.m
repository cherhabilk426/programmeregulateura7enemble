function [idsref,didsref]=Floue_phr(idsref,ephr,vephr);

a1=0.25;a2=0.50;a3=1;

% fonction mu (NG,NM,NP,ZE,PP,PM,PG) pour l'erreur

mueNG=max(min(1,-(ephr+a2)/(a3-a2)),0);
mueNM=max(min((ephr+a3)/(a3-a2),-(ephr+a1)/(a2-a1)),0);
mueNP=max(min((ephr+a2)/(a2-a1),-ephr/a1),0);
mueZE=max(min(1+ephr/a1,1-ephr/a1),0);
muePP=max(min(ephr/a1,(a2-ephr)/(a2-a1)),0);
muePM=max(min((ephr-a1)/(a2-a1),(a3-ephr)/(a3-a2)),0);
muePG=max(min((ephr-a2)/(a3-a2),1),0);

%fonction mu (NG,NM,NP,ZE,PP,PM,PG) pour la dérivée

muveNG=max(min(1,-(vephr+a2)/(a3-a2)),0);
muveNM=max(min((vephr+a3)/(a3-a2),-(vephr+a1)/(a2-a1)),0);
muveNP=max(min((vephr+a2)/(a2-a1),-vephr/a1),0);
muveZE=max(min(1+vephr/a1,1-vephr/a1),0);
muvePP=max(min(vephr/a1,(a2-vephr)/(a2-a1)),0);
muvePM=max(min((vephr-a1)/(a2-a1),(a3-vephr)/(a3-a2)),0);
muvePG=max(min((vephr-a2)/(a3-a2),1),0);

% defuzzification : max-min

musNG=max([min(mueNG,muveNG);min(mueNG,muveNM);min(mueNG,muveNP);
           min(mueNG,muveZE);min(mueNM,muveNG);min(mueNM,muveNM);
           min(mueNM,muveNP);min(mueNP,muveNG);min(mueNP,muveNM);min(mueZE,muveNG)]);
musNM=max([min(mueNG,muvePP);min(mueNM,muveZE);min(mueNP,muveNP);
           min(mueZE,muveNM);min(muePP,muveNG)]);
musNP=max([min(mueNG,muvePM);min(mueNM,muvePP);min(mueNP,muveZE);
           min(mueZE,muveNP);min(muePP,muveNM);min(muePM,muveNG)]);
musZE=max([min(mueNG,muvePG);min(mueNM,muvePM);min(mueNP,muvePP);
           min(mueZE,muveZE);min(muePP,muveNP);min(muePM,muveNM);min(muePG,muveNG)]);
musPP=max([min(mueNM,muvePG);min(mueNP,muvePM);min(mueZE,muvePP);
           min(muePP,muveZE);min(muePM,muveNP);min(muePG,muveNM)]);
musPM=max([min(mueNP,muvePG);min(mueZE,muvePM);min(muePP,muvePP);
           min(muePM,muveZE);min(muePG,muveNP)]);
musPG=max([min(mueZE,muvePG);min(muePP,muvePG);min(muePP,muvePM);
           min(muePM,muvePG);min(muePM,muvePM);min(muePM,muvePP);
           min(muePG,muvePG);min(muePG,muvePM);min(muePG,muvePP);min(muePG,muveZE)]);

som_mus=musNG+musNM+musNP+musZE+musPP+musPM+musPG;

%Id=5;Gud=0.95;

Id=1;

s_PG=1*Id;
s_PM=0.5*Id;
s_PP=0.25*Id;
s_ZE=0*Id;
s_NP=-0.25*Id;
s_NM=-0.5*Id;
s_NG=-1*Id;

didsref=(s_PG*musPG+s_PM*musPM+s_PP*musPP+s_ZE*musZE+s_NP*musNP+s_NM*musNM+s_NG*musNG)/som_mus;

Gud=4.5;
idsref=idsref+Gud*didsref;

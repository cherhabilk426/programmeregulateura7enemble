function xprime=derives(x,vds,vqs,ws,cr);

rs=4.85;rr=3.805;ls=0.274;lr=0.274;lm=0.258;
j=0.031;p=2;f=0.008;

tr=lr/rr;
ts=ls/rs;
s=1-lm*lm/(lr*ls);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phds=x(1);
phqs=x(2);
ids=x(3);
iqs=x(4);
wm=x(5);
thetas=x(6);

dphds=ws*phqs-rs*ids+vds;
dphqs=-ws*phds-rs*iqs+vqs;
ids1=phds/(s*tr*ls)+wm*phqs/(s*ls)-(tr+ts)*ids/(s*ts*tr);
dids=ids1+(ws-wm)*iqs+vds/(s*ls);
iqs1=-wm*phds/(s*ls)+phqs/(s*tr*ls)-(ws-wm)*ids;
diqs=iqs1-(tr+ts)*iqs/(s*ts*tr)+vqs/(s*ls);
dwm=p*p*(iqs*phds-ids*phqs)/j-p*cr*sign(wm)/j-f*wm/j;
dthetas=ws;

xprime=[dphds;dphqs;dids;diqs;dwm;dthetas];


function Tensionref=ondhyst(ias,ibs,ics,ifs,ies,iasref,ibsref,icsref,ifsref,iesref,f1,f2,f3,f4,f5)

pas=0.1;

e1=ias-iasref;
e2=ibs-ibsref;
e3=ics-icsref;
e4=ifs-ifsref;
e5=ies-iesref;

f1_old=f1;

    if e1<-pas
        
      f1=1;
      
   elseif e1>pas
       
      f1=0;
      
   else 
       
      f1=f1_old;
      
   end
   
   f2_old=f2;
   
   if e2<-pas
       
      f2=1;
      
   elseif e2>pas
       
      f2=0;
      
   else 
       
      f2=f2_old;
      
   end
   
   f3_old=f3;
   
   if e3<-pas
       
      f3=1;
      
   elseif e3>pas
       
      f3=0;
      
   else 
       
      f3=f3_old;
      
   end
 
 f4_old=f4;
   
   if e4<-pas
       
      f4=1;
      
   elseif e4>pas
       
      f4=0;
      
   else 
       
      f4=f4_old;
      
   end  
   
 f5_old=f5;
   
   if e5<-pas
       
      f5=1;
      
   elseif e5>pas
       
      f5=0;
      
   else 
       
      f5=f5_old;
      
   end  
   
   
   
vs=514.6;   
va=vs*(2*f1-f2-f3-f4-f5)/5;
vb=vs*(-f1+2*f2-f3-f4-f5)/5;
vc=vs*(-f1-f2+2*f3-f4-f5)/5;
vf=vs*(-f1-f2-f3+2*f4-f5)/5;
ve=vs*(-f1-f2-f3-f4+2*f5)/5;

Tensionref=[va;vb;vc;vf;ve;f1;f2;f3;f4;f5];
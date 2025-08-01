function Bdpt=Boundnode(xscal,yscal,sss)
d=xscal;
d2=yscal;

 Bdpt1=0:d/sss:d;
 Bdpt2=d/sss:d/sss:d-d/sss;
 Bdpt=[0*ones(sss+1,1),d2*Bdpt1';
     d*ones(sss+1,1),d2*Bdpt1';
     Bdpt2',0*ones(sss-1,1);
     Bdpt2',d2*ones(sss-1,1)];
 %plot(Bdpt(:,1),Bdpt(:,2),'.')

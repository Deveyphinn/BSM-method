function [Cenpt NNN]=Collnode(xscal,yscal,s)
xscal=1;
d2=yscal;
x1=0+d/(3*s):d/(3*s):d-d/(3*s);
x2=0+d2/s:d2/s:d2-d2/s;
[X,Y]=meshgrid(x1,x2);

Cenpt=[X(:),Y(:)];
NNN=length(Cenpt);

%plot(Cenpt(:,1),Cenpt(:,2),'.')

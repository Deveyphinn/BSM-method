function [Cenpt NNN]=Centernode(xscal,yscal,s)

d=xscal;
d2=yscal;
x1=0+d/(5*s):d/(5*s):d-d/(5*s);
x2=0+d2/(5*s):d2/(5*s):d2-d2/(5*s);
[X,Y]=meshgrid(x1,x2);

Cenpt=[X(:),Y(:)];
NNN=length(Cenpt);

%plot(Cenpt(:,1),Cenpt(:,2),'.')

%%Code for The discrete-space problem: The BSM method
clear all
clc
warning off; 
format long
syms t
syms x
syms y
%%  
coll=28;%  Number of Coll node,
cen=2;%Number of Center node , OrdeOfM M
bd=36;%Number of Bound node K_1
KP=100; %Number of TBS
[k1,k2]=kkk(KP);%kkk(K)


tic
dt=1/100; 
T=1;
tt=0:dt:T;
tloop=length(tt);
RBFkind=1;%1=x^delta_m, 2=Mq
sp=1.02; %shape parameter
alpha=5; %Scalling parameter

%%
xscal = 1; yscal = 0.5; 

%% define the boundary nodes
Bdpt=Boundnode(xscal,yscal,bd);
%% define the Center nodes
[Cenpt,NumCNode]=Centernode(xscal,yscal,cen);
%%  define the collocation nodes
[Collpt,OrdeOfM]=Collnode(xscal,yscal,coll);
%% Exact function
%The multiplication of the density and the heat capacity
rhocp=@(x,y) 8172*exp(x+y).*520.*exp(x+y);
%% the thermal conductivity i
k11p=@(x,y,t) exp(x+2*y);
k12p=@(x,y,t) 0.5*exp(x+2*y);
k21p=@(x,y,t) 0.5*exp(x+2*y);
k22p=@(x,y,t) 0.8*exp(x+2*y);

k11pdx=@(x,y,t) exp(x+2*y);
k12pdx=@(x,y,t) 0.5*exp(x+2*y);
k21pdy=@(x,y,t) 1*exp(x+2*y);
k22pdy=@(x,y,t) 1.6*exp(x+2*y);
%% the temperature distribution function
Uexact=@(x,y,t)    (2-0.0019*t).*sin(0.5*t).*exp(x+2*y)+50*x+50*y;
dUdxexact=@(x,y,t) (2-0.0019*t).*sin(0.5*t).*exp(x+2*y)+50;
dUdxxexact=@(x,y,t) (2-0.0019*t).*sin(0.5*t).*exp(x+2*y);
dUdyexact=@(x,y,t)2*(2-0.0019*t).*sin(0.5*t).*exp(x+2*y)+50;
dUdyyexact=@(x,y,t)  4*(2-0.0019*t).*sin(0.5*t).*exp(x+2*y);
dUdxyexact=@(x,y,t) 2*(2-0.0019*t).*sin(0.5*t).*exp(x+2*y);
dUdtexact=@(x,y,t) (-0.0019.*sin(0.5*t)+0.5*(2-0.0019*t).*cos(0.5*t)).*exp(x+2*y);
%% the internal heat source
Gfun=@(x,y,t)8172*520*exp(x+y).*exp(x+y).*((-0.0019.*sin(0.5*t)+...
    0.5*(2-0.0019*t).*cos(0.5*t)).*exp(x+2*y))-...
    ( exp(x+2*y).*((2-0.0019*t).*sin(0.5*t).*exp(x+2*y))+...
    0.8*exp(x+2*y).*(4*(2-0.0019*t).*sin(0.5*t).*exp(x+2*y))+...
    (0.5*exp(x+2*y)+0.5*exp(x+2*y)).*(2*(2-0.0019*t).*sin(0.5*t).*exp(x+2*y))+...
    (exp(x+2*y)+1*exp(x+2*y)).*((2-0.0019*t).*sin(0.5*t).*exp(x+2*y)+50)+...
    (0.5*exp(x+2*y)+1.6*exp(x+2*y)).*(2*(2-0.0019*t).*sin(0.5*t).*exp(x+2*y)+50) );


U0=Uexact(Collpt(:,1),Collpt(:,2),0);
dU0dx=dUdxexact(Collpt(:,1),Collpt(:,2),0);
dU0dy=dUdyexact(Collpt(:,1),Collpt(:,2),0);
dU0dxx=dUdxxexact(Collpt(:,1),Collpt(:,2),0);
dU0dyy=dUdyyexact(Collpt(:,1),Collpt(:,2),0);
dU0dxy=dUdxyexact(Collpt(:,1),Collpt(:,2),0);

Gfunn=Gfun(Collpt(:,1),Collpt(:,2),dt);

UN(:,1)=U0;
UNdx(:,1)=dU0dx;
UNdxx(:,1)=dU0dxx;
UNdy(:,1)=dU0dy;
UNdyy(:,1)=dU0dyy;
UNdxy(:,1)=dU0dxy;

%% Defenition of RBF 
if RBFkind==1
phi=@(x,y,xx,yy,m)   ((x-xx).^2+(y-yy).^2+m.^2).^0.5;
dphidx=@(x,y,xx,yy,m) (x-xx)./(((x-xx).^2+(y-yy).^2+m.^2).^0.5);
dphidy=@(x,y,xx,yy,m) (y-yy)./(((x-xx).^2+(y-yy).^2+m.^2).^0.5);
dphidxx=@(x,y,xx,yy,m) 1./(((x-xx).^2+(y-yy).^2+m.^2).^0.5)-(x-xx).^2./(((x-xx).^2+(y-yy).^2+m.^2).^1.5);
dphidyy=@(x,y,xx,yy,m) 1./(((x-xx).^2+(y-yy).^2+m.^2).^0.5)-(y-yy).^2./(((x-xx).^2+(y-yy).^2+m.^2).^1.5);
dphidxy=@(x,y,xx,yy,m) -((2*y - 2*yy).*(x - xx))./(2*((x - xx).^2 + (y - yy).^2 + m^2).^(3/2));    
else
%% Gaussian
phi=@(x,y,xx,yy,m)((x-xx).^2+(y-yy).^2).^6.5;
dphidx=@(x,y,xx,yy,m) 13*(x-xx).*((x-xx).^2+(y-yy).^2).^5.5;
dphidy=@(x,y,xx,yy,m)13*(y-yy).*((x-xx).^2+(y-yy).^2).^5.5;
dphidxx=@(x,y,xx,yy,m)13*((x-xx).^2+(y-yy).^2).^5.5+143*((x-xx).^2).*((x-xx).^2+(y-yy).^2).^4.5;
dphidyy=@(x,y,xx,yy,m)13*((x-xx).^2+(y-yy).^2).^5.5+143*((y-yy).^2).*((x-xx).^2+(y-yy).^2).^4.5;
dphidxy=@(x,y,xx,yy,m) (11*(13*x - 13*xx).*(2*y - 2*yy).*((x - xx).^2 + (y - yy).^2).^(9/2))./2;  
end

%% Defenition of TBF 
Ttheta=@(j,p,x,y,aa) sin(j.*pi.*(x+aa)./2./aa).*sin(p.*pi.*(y+aa)./2./aa);
Ttheta1=@(j,p,x,y,aa) cos(j.*pi.*(x+aa)./2./aa).*sin(p.*pi.*(y+aa)./2./aa).*j.*pi./2./aa;
Ttheta2=@(j,p,x,y,aa) sin(j.*pi.*(x+aa)./2./aa).*cos(p.*pi.*(y+aa)./2./aa).*p.*pi./2./aa;
Ttheta11=@(j,p,x,y,aa) -sin(j.*pi.*(x+aa)./2./aa).*sin(p.*pi.*(y+aa)./2./aa).*(j.*pi./2./aa).^2;
Ttheta22=@(j,p,x,y,aa) -sin(j.*pi.*(x+aa)./2./aa).*sin(p.*pi.*(y+aa)./2./aa).*(p.*pi./2./aa).^2;
Ttheta12=@(j,p,x,y,aa) (j.*p*pi^2.*cos((pi*j.*(aa + x))./(2*aa)).*cos((pi*p.*(aa + y))./(2*aa)))./(4*aa^2);  

theta=Ttheta(k1,k2,Bdpt(:,1),Bdpt(:,2),alpha);
Theta=Ttheta(k1,k2,Collpt(:,1),Collpt(:,2),alpha);
Theta1=Ttheta1(k1,k2,Collpt(:,1),Collpt(:,2),alpha);
Theta2=Ttheta2(k1,k2,Collpt(:,1),Collpt(:,2),alpha);
Theta11=Ttheta11(k1,k2,Collpt(:,1),Collpt(:,2),alpha);
Theta22=Ttheta22(k1,k2,Collpt(:,1),Collpt(:,2),alpha);
Theta12=Ttheta12(k1,k2,Collpt(:,1),Collpt(:,2),alpha);

%% finiding unknown cooefecient in Eq. (24)
qm=theta\(-phi(Bdpt(:,1),Bdpt(:,2),Cenpt(:,1)',Cenpt(:,2)',sp)); 
CollPhi=phi(Collpt(:,1),Collpt(:,2),Cenpt(:,1)',Cenpt(:,2)',sp)+Theta*qm;
ColldPhidx=dphidx(Collpt(:,1),Collpt(:,2),Cenpt(:,1)',Cenpt(:,2)',sp)+Theta1*qm;
ColldPhidy=dphidy(Collpt(:,1),Collpt(:,2),Cenpt(:,1)',Cenpt(:,2)',sp)+Theta2*qm;
ColldPhidxx=dphidxx(Collpt(:,1),Collpt(:,2),Cenpt(:,1)',Cenpt(:,2)',sp)+Theta11*qm;
ColldPhidyy=dphidyy(Collpt(:,1),Collpt(:,2),Cenpt(:,1)',Cenpt(:,2)',sp)+Theta22*qm;
ColldPhidxy=dphidxy(Collpt(:,1),Collpt(:,2),Cenpt(:,1)',Cenpt(:,2)',sp)+Theta12*qm;

rhocpp=rhocp(Collpt(:,1),Collpt(:,2));


k11pp=k11p(Collpt(:,1),Collpt(:,2));
k12pp=k12p(Collpt(:,1),Collpt(:,2));
k21pp=k21p(Collpt(:,1),Collpt(:,2));
k22pp=k22p(Collpt(:,1),Collpt(:,2));
k11pdxp=k11pdx(Collpt(:,1),Collpt(:,2));
k12pdxp=k12pdx(Collpt(:,1),Collpt(:,2));
k21pdyp=k21pdy(Collpt(:,1),Collpt(:,2));
k22pdyp=k22pdy(Collpt(:,1),Collpt(:,2));


for i=1:tloop-1
    i;
    time=dt*i;
  %%finiding unknown cooefecient in Eq. (24)
Pm(:,i)=theta\(Uexact(Bdpt(:,1),Bdpt(:,2),time));
Up=Theta*Pm(:,i); dUpdx=Theta1*Pm(:,i); dUpdy=Theta2*Pm(:,i); 
dUpdxx=Theta11*Pm(:,i); dUpdyy=Theta22*Pm(:,i); dUpdxy=Theta12*Pm(:,i);


UN(:,i)=U0;
 UNdx(:,i)=dU0dx;
UNdxx(:,i)=dU0dxx;
UNdy(:,i)=dU0dy;
UNdyy(:,i)=dU0dyy;
UNdxy(:,i)=dU0dxy;

Gfunn=Gfun(Collpt(:,1),Collpt(:,2),time);
%% the right side of Eq. (23)
f1=Gfunn+rhocpp*(1/(2*dt)).*UN(:,i-1);

LUp=rhocpp*(1/(2*dt)).*Up-...
    ( k11pp.*dUpdxx+k22pp.*dUpdyy+(k12pp+k21pp).*dUpdxy+...
(k11pdxp+k21pdyp).*dUpdx+(k12pdxp+k22pdyp).*dUpdy);

F1=(f1-LUp);
  %%finiding unknown cooefecient in Eq. (25)
q1(:,i)=(rhocpp.*(1/(2*dt)).*CollPhi-...
    ( k11pp.*ColldPhidxx+k22pp.*ColldPhidyy+(k12pp+k21pp).*ColldPhidxy+...
(k11pdxp+k21pdyp).*ColldPhidx+(k12pdxp+k22pdyp).*ColldPhidy))\F1;

%% deassign
U0=     Up    +CollPhi*q1(:,i);
dU0dx=  dUpdx +ColldPhidx*q1(:,i);
dU0dy=  dUpdy +ColldPhidy*q1(:,i);

dU0dxx= dUpdxx+ColldPhidxx*q1(:,i);
dU0dyy= dUpdyy+ColldPhidyy*q1(:,i);
dU0dxy= dUpdxy+ColldPhidxy*q1(:,i);




%% Error
ErrorUmae=max(abs(Uexact(Collpt(:,1),Collpt(:,2),time)-U0));
 Errortime(i)=ErrorUmae;
 
 
end

max(Errortime)






function [K1 K2]=kkk(s)

k1=1;
for i=2:20
    qq=1:i;
    k1=[k1;qq'];
    
end
k1=k1(1:s);k2=1;
for i=2:20
    qq=1:i;
    pp=fliplr(qq);
    k2=[k2;pp'];
   
end
 k2=k2(1:s);
 K1=k1';
 K2=k2';
 
 
 
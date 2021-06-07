function [yy]=rangleconv(nb,nx,xx);
%nb=滤波器长度，nx数据长度，xx输入数据
nv=nx+nb-1;%nv输出数据长度；
for i=1:nv;
    bb(i)=0;
end;
 bb(1)=xx(1);
 for i=2:nx;
     bb(i)=bb(i-1)+xx(i);%B(Z)=X(Z)/(1-Z)
 end;
 for i=nx+1:nv;
     bb(i)=bb(i-1);
 end;
 for i=1:nb;
     yy(i)=bb(i);
 end;
 for i=nb+1:nv;
     yy(i)=bb(i)-bb(i-nb);%YZ=B(Z)*(1-Z^k)
 end;
 for i=1:nv;
     yy(i)=yy(i)/nb;
 end
 
     
function [ vv ] = triangle(nr,m1,n12,uu)
%    uu=ones(1,10);
%    nr=3;
%    m1=1;
%    n12=10;
    yy=rangleconv(nr,n12,uu);
   np=nr+n12-1;
   d=rangleconv(nr,np,yy);


   for i= 1:n12;
                        tt(i) = d(i+nr-1);
   end;
   for i= 1:nr-1  ;                               
                        tt(i) = tt(i) + d(nr-i);
   end;
   for i= 1:nr-1;                                 
                        tt(n12-i+1) = tt(n12-i+1) + d(n12+(nr-1)+i);
   end
  for i=1:n12;
      vv(i)=tt(i);
  end
  function [simi] = Localsimilarity(G,E,N1,N2,mu1,mu2,itermax,sol)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function used for obtaining the predictability of multiple by computing local similarity
% Code by:  Junming Zhang  Deli Wang  Bin Hu
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G:full wavefield data spectrum
% E:predicted multiple spcetrum
% N1:  smooth length in 1st dimensions
% N2  smooth length in 2st dimensions
% mu1:parameters controlling the physical dimensionality
% mu2:parameters controlling the physical dimensionality
% sol ls or sr
  [nt,nh] = size(E);
  scaleE=max(abs(E(:)));
  scaleG=max(abs(G(:)));
  E=E./scaleE;
  G=G./scaleG;
  flt=zeros(2*N1-1,1);
  flt(1:N1)=1:N1;
  flt(N1:2*N1-1)=flt(N1:-1:1);
  lap1=zeros(nt);
  window=nt;
  for ii=1:nt
      if ((ii-N1)<0)
          lap1(1:N1-1+ii,ii)=flt(N1+1-ii:end);
      else if((ii+N1-1)>nt)
              lap1(ii-N1+1:end,ii)=flt(1:nt-ii+N1);
          else 
              lap1(ii-N1+1:ii+N1-1,ii)=flt;
          end
      end
      lap1(:,ii)=lap1(:,ii)./sum(lap1(:,ii));  
  end
  lap1=sparse(lap1);
  lap1t=lap1';

  flt=zeros(2*N2-1,1);
  flt(1:N2)=1:N2;
  flt(N2:2*N2-1)=flt(N2:-1:1);
  nx=nh;
  lap2=zeros(nx);
  for ii=1:nx
      if ((ii-N2)<0)
          lap2(1:N2-1+ii,ii)=flt(N2+1-ii:end);
      else if((ii+N2-1)>nx)
              lap2(ii-N2+1:end,ii)=flt(1:nx-ii+N2);
          else 
              lap2(ii-N2+1:ii+N2-1,ii)=flt;
          end
      end
      lap2(:,ii)=lap2(:,ii)./sum(lap2(:,ii));
  end
    lap2(isnan(lap2)==1)=0;
  lap2=sparse(lap2);
  lap2t=lap2';
  V=zeros(nt,nh);
  
  method=0;
  gap=floor(nh/10);
  for ih=1:nh;
      if mod(ih,gap)==0
          disp(['The completion rate:',num2str(ih/gap*10),'%;']);
      end       
           a=E(:,ih);
           b=G(:,ih); 
           A=sparse(1:nt,1:nt,a,nt,nt);
           B=sparse(1:nt,1:nt,b,nt,nt);
           zeo1=sparse(1:nt,1:nt,a.*a,nt,nt);
           zeo2=sparse(1:nt,1:nt,b.*b,nt,nt);
           eig1=eig(zeo1);
           eig2=eig(zeo2);
           [var index1]=max(gradient(eig1));
         
            if eig1(index1)==0;
               index1=index1+1;
           end
           [var index2]=max(gradient(eig2));
           if eig2(index2)==0;
               index2=index2+1;
           end
           if isequal(sol,'sr'); 
               d1=mu1*eig1(index1);
               d2=mu2*eig2(index2);
               c1l=d1*speye(nt)+lap1t*(zeo1-d1*speye(nt))*lap1;
               c2l=d2*speye(nt)+lap1t*(zeo2-d2*speye(nt))*lap1;
           else isequal(sol,'ls');
               d1=mu1*eig1(index1);
               d2=mu2*eig2(index2);               
               c2l=zeo1+eig1(index1)*speye(window)*nt;
               c1l=zeo2+eig2(index2)*speye(window)*nt;
           end
           c1r=(lap1t*B)*a;
           %BTa
           c2r=(lap1t*A)*b; 
           if ih==1
           tic
           EE=speye(nt);
           X0=sparse(1:nt,1:nt,1./diag(c1l),nt,nt);    
           for iter=1:itermax
               Xn=X0*(2*EE-c1l*X0); 
               X0=Xn; 
           end
           c1=lap1*(Xn*c1r);
           timeiter=toc;
           tic
           c1=lap1*(((c1l)\eye(nt))*c1r);
           timeinverse=toc;
           if (timeinverse<timeiter); method=1;end
           end
           if (method==1)
               c1=lap1*(((c1l)\eye(nt))*c1r);
               c2=lap1*(((c2l)\eye(nt))*c2r);
           else
               EE=speye(nt);
               X0=sparse(1:nt,1:nt,1./diag(c1l),nt,nt);   
               for iter=1:itermax
                   Xn=X0*(2*EE-c1l*X0); 
                   X0=Xn; 
               end
                c1=lap1*(Xn*c1r);
               X0=sparse(1:nt,1:nt,1./diag(c2l),nt,nt);   
               for iter=1:itermax
                   Xn=X0*(2*EE-c2l*X0); 
                   X0=Xn; 
               end
               c2=lap1*(Xn*c2r);
           end
           V1(:,ih)=c1;
           V2(:,ih)=c2;
  end; 
  
  for ii=1:itermax  
  V2=V2*lap2;
  V1=V1*lap2;
  end


V1(isnan(V1)==1)=0;
V2(isnan(V2)==1)=0;
simi=V1.*V2;
  aa=find(simi<0);
  simi(aa)=0;
  simi=simi./max(simi(:));
% end


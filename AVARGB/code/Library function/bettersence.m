 function [R,G,B] = bettersence(R,G,B,fR11,pickpoint)
% pickpoint=fR11;
is=1;
      sence=ones(1,3);
      for it1=1:length(fR11)
              if fR11(it1,4)<50&&fR11(it1,5)<160&&fR11(it1,6)<160
                  sence(is,1:2)=fR11(it1,1:2);
                  sence(is,3)=it1;
                  is=is+1;
              end
      end
                 [ls1,ws1]=size(sence);
                 [ls2,ws2]=size(pickpoint);
                 for is1=1:ls1
                     for is2=1:ls2
                         if sence(is1,1)==pickpoint(is2,1)&&sence(is1,2)==pickpoint(is2,2)
                             sence(is1,:)=0;
                         end
                     end
                 end
                 [l0,w0]=find(sence(:,1)==0);
                 sence(l0,:)=[];
                 fR11(sence(:,3),4:6)=255;
                 for ifr=1:length(fR11)   
                 R(fR11(ifr,1),fR11(ifr,2))=fR11(ifr,4);
                 G(fR11(ifr,1),fR11(ifr,2))=fR11(ifr,5);
                 B(fR11(ifr,1),fR11(ifr,2))=fR11(ifr,6);
                 end
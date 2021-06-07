 function [shape,pl2,pickpoint] = Pick(v1,RGB,fR11,dt,vmin,dv,limit,pl,pv,pa)
 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function used for picking primary energy groups
% Code by:  Junming Zhang, Bin Hu, Deli Wang and Xiangbo Gong 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v1: The velocity spectrum of full wavefield data
% RGB: RGB space
% fR11: The RGB value of energy groups
% dt: time sample 
% vmin: Minimum scan velocity
% dv: Scan velocity step
% limit: Y,N,Whether the selection limit of each feature is needed to constrain the primary pickup
% pl:Upper limit of predictability pickup
% pv:Upper limit of velocity pickup
% pa:Upper limit of amplitude pickup
pickpoint=fR11;
R=RGB(:,:,1);
G=RGB(:,:,2);
B=RGB(:,:,3);
 pickpoint(:,3)=pickpoint(:,4).*5+pickpoint(:,5).*2+pickpoint(:,6)*2;
 pickpoint=sortrows(pickpoint,3);
[Lpickpoint,Wpickpoint]=find(pickpoint(:,1)==0);
pickpoint(Lpickpoint,:)=[];
cpoint=pickpoint(:,3);
cpoint1=gradient(cpoint,1);
[Lcpoint1,Wcpoint1]=size(cpoint1);
cpoint2=cpoint1(1:floor(Lcpoint1/2));
maxnumber=max(cpoint);
[aa,Num]=max(cpoint2);
choosenumber=maxnumber-cpoint(Num);
tpickpoint=pickpoint;
[Ltpickpoint,Wpickpoint]=size(tpickpoint);
for itpickpoint=1:Ltpickpoint
    if maxnumber-tpickpoint(itpickpoint,3)<choosenumber
        tpickpoint(itpickpoint,:)=0;
    end
end
[Ltpickpoint,Wpickpoint]=find(tpickpoint(:,1)==0);
tpickpoint(Ltpickpoint,:)=[];
tpickpoint=sortrows(tpickpoint,1);
[Ltpickpoint,Wtpickpoint]=size(tpickpoint);

%% useing selection limit of each feature to constrain the primary pickup 

if strcmp(limit,'Y')
for ip1=1:Ltpickpoint
    if tpickpoint(ip1,4)>=pl
        fillpoint(ip1,:)=tpickpoint(ip1,:);
        tpickpoint(ip1,:)=0;
        
    elseif tpickpoint(ip1,5)>=pv
         fillpoint(ip1,:)=tpickpoint(ip1,:);
         tpickpoint(ip1,:)=0;
     elseif tpickpoint(ip1,6)>=pa;
          fillpoint(ip1,:)=tpickpoint(ip1,:);
          tpickpoint(ip1,:)=0;
    elseif tpickpoint(ip1,:)<0.3/dt
        tpickpoint(ip1,:)=0;
    end
   
end
[Ltpickpoint,Wtpickpoint]=find(tpickpoint(:,1)==0);
tpickpoint(Ltpickpoint,:)=[];
end

% % Remove energy groups that do not meet the pickup rules

[CLtpickpoint,CWtpickpoint]=size(tpickpoint);
if CLtpickpoint>0
ttpickpoint=tpickpoint(:,1:2);
ffpickpoint=ttpickpoint;
 for tpickpoint1=1:CLtpickpoint
           TWmin=tpickpoint(tpickpoint1,1)-0.15/dt;
           TWmax=tpickpoint(tpickpoint1,1)+0.15/dt;
      for tpickpoint2=1:CLtpickpoint        
       if tpickpoint(tpickpoint2,1)>TWmin&&tpickpoint(tpickpoint2,1)<TWmax&&tpickpoint(tpickpoint2,1)~=tpickpoint(tpickpoint1,1)
           if v1(tpickpoint(tpickpoint2,1),tpickpoint(tpickpoint2,2))<v1(tpickpoint(tpickpoint1,1),tpickpoint(tpickpoint1,2))
               ffpickpoint(tpickpoint2,:)=0;
           else ffpickpoint(tpickpoint1,:)=0;
           end
       end
      end
 end
 
 tpickpoint=ffpickpoint;
    [Lpickpoint,Wpickpoint]=find(tpickpoint(:,1)==0);
tpickpoint(Lpickpoint,:)=[];

end
[CLtpickpoint,CWtpickpoint]=size(tpickpoint);
if CLtpickpoint>3
 pickpoint=tpickpoint;
else
lP=4-CLtpickpoint;
addPick=fillpoint(1:lP,1:2);
pickpoint=tpickpoint(:,1:2);
pickpoint=[pickpoint;addPick];
end

%%Fitting the picked primary  to the curve

pickpoint=sortrows(pickpoint,1);

 [tpickpoint,vpickpoint]=find(pickpoint(:,1)>0);
 for i=1:length(tpickpoint)
     fpickpoint(i,:)=pickpoint(tpickpoint(i,1),:);
 end
 pickpoint=fpickpoint(:,1:2);
 fp=interp1(pickpoint(:,1),pickpoint(:,2),min(pickpoint(:,1)):max(pickpoint(:,1)),'cubic');
 Fx2=gradient(fp,1);
 Fx2=Fx2';
 Fx2=[zeros(pickpoint(1,1)-1,1);Fx2];
 for ifp=1:length(pickpoint)
    if Fx2(pickpoint(ifp,1),1)<0
        break
    end
end
if ifp<3
    ifp=3;
end
fpick1=pickpoint(1:ifp,:);
fpick2=pickpoint(ifp:length(pickpoint),:);

  pg31=polyfit(fpick1(:,1),fpick1(:,2),2);
   pg32=polyfit(fpick2(:,1),fpick2(:,2),2);
  tfinal=length(R);
   vfinal=tfinal^2*pg32(:,1)+tfinal*pg32(:,2)+pg32(:,3);
  tstart=1;
      vstart=(1430-vmin)/dv;
 pickstart=[tstart,vstart];
 pickfinal=[tfinal,vfinal];
 pickpoint=[pickstart;pickpoint(:,1:2);pickfinal];
 pl1=interp1(pickpoint(:,1),pickpoint(:,2),min(pickpoint(:,1)):max(pickpoint(:,1)),'cubic');
  pl2=triangle(10,1,length(pl1),pl1);

     [nr,nc]=size(R);
     for ipickpoint=1:length(pickpoint)
         for itpickpoint=1:length(ttpickpoint)
             if ttpickpoint(itpickpoint,:)==pickpoint(ipickpoint,:)
                 ttpickpoint(itpickpoint,:)=0;
             end
         end
     end
     [Ltpickpoint,Wtpickpoint]=find(ttpickpoint(:,1)==0);
      ttpickpoint(Ltpickpoint,:)=[];
    [Lttpickpoint,Wttpickpoint]=size(ttpickpoint);
      for itpickpoint2=1:Lttpickpoint
         for ifR11=1:length(fR11)
             if fR11(ifR11,1)==ttpickpoint(itpickpoint2,1)&&fR11(ifR11,2)==ttpickpoint(itpickpoint2,2)
                 fR11(ifR11,:)=0;
             end
         end
     end
     [LfR11,WfR11]=find(fR11(:,1)==0);
      fR11(LfR11,:)=[];
      
     %%Highlight the color of each energy group
tya=10;
tyb=1;
tyc=floor((tya^2-tyb^2)^0.5);
RR=R;
GG=G;
BB=B;
for rc=1:length(fR11)
cc=[fR11(rc,1),fR11(rc,2)];
cc1=[cc(:,1)-tyc,cc(:,2)];
cc2=[cc(:,1)+tyc,cc(:,2)];
for i=1:nr
for j=1:nc
t=[i j];
ls1=norm(cc1-t);
ls2=norm(cc2-t);
ls=ls1+ls2;
   if ls<2*tya
RR(i,j)=R(fR11(rc,1),fR11(rc,2));
GG(i,j)=G(fR11(rc,1),fR11(rc,2));
BB(i,j)=B(fR11(rc,1),fR11(rc,2));
shape=cat(3,RR,GG,BB); 
   end
end
end
end

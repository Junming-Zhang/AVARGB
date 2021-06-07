                          function [RGB,gg] = RGBmapping(v1,simi,N1,dt,vmin,vmax,nv,maxpick,border)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function used for building RGB space
% Code by:  Junming Zhang, Bin Hu, Deli Wang and Xiangbo Gong 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v1: The velocity spectrum of full wavefield data
% simi: Similarity  spectrum of full wavefield data and predicted multiple
% N1: The length of the window for calculating average predictability
% dt: time sample 
% vmin: Minimum scan velocity
% vmax: maxium scan velocity
% nv: scan number
% maxpick:The length of local maxima window
% border: Y,N Is it necessary to exclude energy groups with  abnormally high and low velocities

dv=(vmax-vmin)/nv;
   tmx=mean(v1(:));
   [nt,nh]=size(v1);
   gg=zeros(nt,6);
   mv1=zeros(size(v1));
   for imx=1:nh
   z=v1(:,imx);
   [pks2,locs2]=findpeaks(z,'MinPeakHeight',tmx);%0.001
   mv1(locs2,imx)=v1(locs2,imx);
   end
ms=medfilt3(v1,[12,5]);
   mx=matrix_local_maximum_improved1(ms,maxpick);
    mx=triangle2(2,3,nh,nt,mx);
   mx1=matrix_local_maximum_improved1(mv1,maxpick);
    mx1=triangle2(2,3,nh,nt,mx1);
   mx2=mx.*mx1;
   mx3=matrix_local_maximum_improved1(mx2,maxpick);

 %%RGB
R=255*ones(size(mv1));
G=255*ones(size(mv1));
B=255*ones(size(mv1));
R=uint8(R);
G=uint8(G);
B=uint8(B);
    [tm,vm]=find(mx3>0);
    P1=[tm,vm];
    P1t=sortrows(P1,2);
    [LP1,WP1]=size(P1t);
    %%exclude energy groups with  abnormally high and low velocities
   if strcmp(border,'Y')
    for iLP1=1:LP1
        if P1t(iLP1,2)==1
            P1t(iLP1,:)=0;
        elseif P1t(iLP1,2)==nv
            P1t(iLP1,:)=0;
        end
    end
   [LTP1,WTP1]=find(P1t(:,1)==0);
   P1t(LTP1,:)=[];
   [nP1,mP1]=size(P1t);
   nboder=floor(0.96*nP1);
   deboder=P1t(nboder,2);
   for inP1=1:nP1
       if P1t(inP1,2)>deboder
           P1t(inP1,:)=0;
       end
   end
    [LTP1,WTP1]=find(P1t(:,1)==0);
   P1t(LTP1,:)=[];
   end

    %%     Convert multiple predictability into R-axis coordinates   
for  ii=1:length(P1t);
             if ((ii-floor(N1/2))<0)
             lwindow=simi(1:P1t(ii,1)+floor(N1/2)-1,P1t(ii,2));
      else if(ii+floor(N1/2)-1)>length(simi)
             lwindow=simi(P1t(ii,1)-floor(N1/2)+1:end,P1t(ii,2));
          else 
             lwindow=simi(P1t(ii,1)-floor(N1/2)+1:P1t(ii,1)+floor(N1/2)-1,P1t(ii,2));
          end
             end      
        local(ii)=mean(lwindow(:));
end


% % % % % % % % % 
  R1=[P1t,local'];
  R1t=sortrows(R1,1);
  for iR1t=2:length(R1t)-1
     if R1t(iR1t-1,1)<0.3/dt;
         R1t(iR1t-1,:)=0;
     end
  end
  tR1t=find(R1t(:,1)>0);
  for itR1t=1:length(tR1t)
      R11t(itR1t,:)=R1t(tR1t(itR1t,1),:);
  end
  R11=sortrows(R11t,3);

r1tt=log10(R11(:,3));
r1tts=triangle(ceil(length(r1tt)/5),1,length(r1tt),r1tt);
Fxr1tts=gradient(r1tts,1);
 
Fxr2tts=gradient(Fxr1tts,1);
Fxr2tts=Fxr2tts';
% save Fxr2tts
thr=zeros(1,1);
for iFx=2:length(Fxr2tts)
    if Fxr2tts(iFx-1,1)*Fxr2tts(iFx,1)<0
        thr(iFx,1)=iFx;
    end
end
thr=deletesim(thr);
if length(thr)>=2;
    thr1=R11(thr(1,1),3);
    thr2=R11(thr(length(thr),1),3);
elseif length(thr)<2
    thr1=R11(thr(1,1),3);
    meanthr=mean(abs(Fxr2tts));
    for imean=2:length(Fxr2tts)-1
    if abs(Fxr2tts(imean-1))<meanthr&&abs(Fxr2tts(imean))>meanthr
        break
    elseif abs(Fxr2tts(imean-1))>meanthr&&abs(Fxr2tts(imean))<meanthr
    break
    end
    end
    thr2=R11(imean,3);
    thr=[thr;imean];
end
 thrmix=R11(thr(1,1)+1:thr(length(thr),1)-1,:);  
 [Xthmix,Ythmix]=size(thrmix);

 mixmin=abs(floor(log10(thrmix(1,3))));
 mixmax=abs(floor(log10(thrmix(Xthmix,3))));
 thr3=floor((abs(mixmax-mixmin)+2)/2);
thrmix1=thr1*10^thr3;
thrmix2=thr2/10^thr3;
Rmix=255/3;
 X1=[thr1,thrmix1];
 Y1=[1,Rmix];
 X2=[thrmix1,thrmix2];
 Y2=[Rmix,2*Rmix];
 X3=[thrmix2,thr2];
 Y3=[2*Rmix,255];
p1= polyfit(X1,Y1,1);
p2= polyfit(X2,Y2,1); 
p3= polyfit(X3,Y3,1);

   for im=1:length(R11)
       if R11(im,1)<=0.4/dt;
           R2=255;
       elseif R11(im,3)<=thr1
           R2=0;
       elseif R11(im,3)>=thr2;
           R2=255;
       else
           if R11(im,3)<thrmix1&&R11(im,3)>thr1

               R2=R11(im,3)*p1(:,1)+p1(:,2);
                 
           elseif R11(im,3)<thrmix2&&R11(im,3)>thrmix1
                R2=R11(im,3)*p2(:,1)+p2(:,2);
           elseif R11(im,3)<thr2&&R11(im,3)>thrmix2
                R2=R11(im,3)*p3(:,1)+p3(:,2);
           end   
       end
       R(R11(im,1),R11(im,2))=R2;
   end  
 
 
%% Convert velocity to G-axis coordinates
ht=floor(nt/5);
[tR12,vR12]=find(R<255);
R12=[tR12,vR12];
tGx1=find(R12(:,1)<ht);
if ~isempty(tGx1)
for itGx1=1:length(tGx1)
TGX1(itGx1,:)=R12(tGx1(itGx1,:),:);
end
else
    tGx1=find(R11(:,1)<ht);
for itGx1=1:length(tGx1)
TGX1(itGx1,:)=R11(tGx1(itGx1,:),:);
end
end
tGx2=find(R12(:,1)>ht&R12(:,1)<2*ht);
if ~isempty(tGx2)
for itGx2=1:length(tGx2)
TGX2(itGx2,:)=R12(tGx2(itGx2,:),:);
end
else
    tGx2=find(R11(:,1)>ht&R11(:,1)<2*ht);
for itGx2=1:length(tGx2)
TGX2(itGx2,:)=R11(tGx2(itGx2,:),:);
end
end
tGx3=find(R12(:,1)>2*ht&R12(:,1)<3*ht);
if length(tGx3)>0
for itGx3=1:length(tGx3)
TGX3(itGx3,:)=R12(tGx3(itGx3,:),:);
end
else
    tGx3=find(R11(:,1)>2*ht&R11(:,1)<3*ht);
for itGx3=1:length(tGx3)
TGX3(itGx3,:)=R11(tGx3(itGx3,:),:);
end
end
tGx4=find(R12(:,1)>3*ht&R12(:,1)<4*ht);
if ~isempty(tGx4)
for itGx4=1:length(tGx4)
TGX4(itGx4,:)=R12(tGx4(itGx4,:),:);
end
else 
    tGx4=find(R11(:,1)>3*ht&R11(:,1)<4*ht);
for itGx4=1:length(tGx4)
TGX4(itGx4,:)=R11(tGx4(itGx4,:),:);
end
end
tGx5=find(R12(:,1)>4*ht&R12(:,1)<5*ht);
if ~isempty(tGx5)
for itGx5=1:length(tGx5)
TGX5(itGx5,:)=R12(tGx5(itGx5,:),:);
end
else
    tGx5=find(R11(:,1)>4*ht&R11(:,1)<5*ht);
for itGx5=1:length(tGx5)
TGX5(itGx5,:)=R11(tGx5(itGx5,:),:);
end
end

TGX1=sortrows(TGX1,2);
[l1,w1]=size(TGX1);
TGX2=sortrows(TGX2,2);
[l2,w2]=size(TGX2);
TGX3=sortrows(TGX3,2);
[l3,w3]=size(TGX3);
TGX4=sortrows(TGX4,2);
[l4,w4]=size(TGX4);
TGX5=sortrows(TGX5,2);
[l5,w5]=size(TGX5);
Gx1=TGX1(1,:);
Gx2=TGX2(l2,:);
Gx3=TGX3(l3,:);
Gx4=TGX4(l4,:);
Gx5=TGX5(l5,:);


gX1=[Gx1(:,1),Gx2(:,1),Gx3(:,1),Gx4(:,1),Gx5(:,1)];
gY1=[Gx1(:,2),Gx2(:,2),Gx3(:,2),Gx4(:,2),Gx5(:,2)];
 
 Gs=(1400-vmin)/dv; 
 yX1=[1,gX1];
 yY1=[Gs,gY1];
 Yq=interp1(yX1,yY1,1:Gx5(:,1),'linear');
 pg2=polyfit(1:length(Yq),Yq,1);
 Gf=nt*pg2(:,1)+pg2(:,2);
gX=[1,gX1,nt];
gY=[Gs,gY1,Gf];
 tR11=R11(:,1);
vR11=R11(:,2);
G1=[tR11,vR11];
yq=interp1(gX,gY,1:nt,'linear');
    yq1=triangle(200,1,nt,yq); 

  for igx=1:length(R11)
    Gx(igx)=vR11(igx)/yq1(tR11(igx));
end

Gx6=[0.5,1];
Gy=[0,255];
G1=[G1,Gx'];
 G1=sortrows(G1,1); 

 pg=polyfit(Gx6,Gy,1);
for lm=1:length(G1);
    if G1(lm,3)>1
        G2=255;
    else if G1(lm,3)<0.5
            G2=0;
        else
    G2=G1(lm,3)*pg(:,1)+pg(:,2);
        end
    end
    G(G1(lm,1),G1(lm,2))=255-G2;
end

%%Convert amplitude to B-axis coordinates

twindow=floor(0.15/dt);
B1=G1(:,1:2);
[bmax,bmin]=amplitude(v1,twindow);
for ib=1:nt
 bmixmin=abs(floor(log10(bmin(ib))));
 bmixmax=abs(floor(log10(bmax(ib)))); 
 bhr=floor((abs(bmixmax-bmixmin))/3);
bmix1(ib)=bmin(ib)*10^bhr;
bmix2(ib)=bmax(ib)/10^bhr;
end
Bmix=255/3;
for ib1=1:length(B1)
 X1=[bmin(B1(ib1,1)),bmix1(B1(ib1,1))];
 Y1=[255,2*Bmix];
 X2=[bmix1(B1(ib1,1)),bmix2(B1(ib1,1))];
 Y2=[2*Bmix,Bmix];
 X3=[bmix2(B1(ib1,1)),bmax(B1(ib1,1))];
 Y3=[Bmix,1];
bpg1= polyfit(X1,Y1,1);
bpg2= polyfit(X2,Y2,1); 
bpg3= polyfit(X3,Y3,1);
if v1(B1(ib1,1),B1(ib1,2))<bmix1(ib1)
    B1(ib1,3)=v1(B1(ib1,1),B1(ib1,2))*bpg1(:,1)+bpg1(:,2);
elseif bmix1(ib1)<=v1(B1(ib1,1),B1(ib1,2))&&v1(B1(ib1,1),B1(ib1,2))<bmix2(ib1)
     B1(ib1,3)=v1(B1(ib1,1),B1(ib1,2))*bpg2(:,1)+bpg2(:,2);
elseif bmix2(ib1)<ms(B1(ib1,1),B1(ib1,2))
     B1(ib1,3)=v1(B1(ib1,1),B1(ib1,2))*bpg3(:,1)+bpg3(:,2);
end
B(B1(ib1,1),B1(ib1,2))=B1(ib1,3);
B1(ib1,4)=v1(B1(ib1,1),B1(ib1,2));
end

%%Remove abnormal points
 ij=1;
for ii=1:length(G1);
    G1(ii,4)=R(G1(ii,1),G1(ii,2));
    G1(ii,5)=G(G1(ii,1),G1(ii,2));
    G1(ii,6)=B(G1(ii,1),G1(ii,2));
    if R(G1(ii,1),G1(ii,2))<150&&G(G1(ii,1),G1(ii,2))<160&&B(G1(ii,1),G1(ii,2))<220;
   T1(ij,1:2)=G1(ii,1:2);
   ij=ij+1;
    end
end
   T1=sortrows(T1,1);  

[LT1,WT1]=size(T1);
T2=T1;
T3=zeros(LT1,1,LT1);
for iT1=1:LT1
    TWmin=T1(iT1,1)-0.2/dt;
    TWmax=T1(iT1,1)+0.2/dt;
    for IT1=1:LT1   
        if T1(IT1,1)>=TWmin&&T1(IT1,1)<=TWmax
            T3(IT1,1,iT1)=ms(T1(IT1,1),T1(IT1,2));
        end
    end
     if ms(T1(iT1,1),T1(iT1,2))==max(T3(:,:,iT1));
         T2(iT1,:)=0;

     end
end
fR1=T1-T2;
[Lf1,Wf1]=find(fR1(:,1)==0);
fR1(Lf1,:)=[];
[LfR1,WfR1]=size(fR1);
tfR1=zeros(1,2);
for ifR1=1:LfR1
    fR1(ifR1,3)=vmin+dv*fR1(ifR1,2);
      fR1(ifR1,4)=0.02*fR1(ifR1,3);
end

T2=[T2;tfR1];
fj2=1;
tfR2=zeros(1,2);
for ifR2=1:LfR1
   if fR1(ifR2,1)>nt-50;
     tfR2(fj2,:)=fR1(ifR2,1:2);
            fj2=fj2+1;   
   end
end

if tfR2(1,1)>0  
    T2=[T2;tfR2];
end

[LT2,WT2]=find(T2(:,1)==0);
T2(LT2,:)=[]; 
   T4=T2;
[LT4,WT4]=size(T4);
for iT2=1:LT4   
        R(T4(iT2,1),T4(iT2,2))=255;
        G(T4(iT2,1),T4(iT2,2))=255;
        B(T4(iT2,1),T4(iT2,2))=255;
end

fR11=G1;
if LT4>0
for iG1=1:length(fR11)
    for iT4=1:LT4
        if fR11(iG1,1)==T4(iT4,1)&&fR11(iG1,2)==T4(iT4,2);
            fR11(iG1,:)=0;
        end
    end
end
[LG1,WG1]=find(fR11(:,1)==0);
fR11(LG1,:)=[];
end
for iG2=1:length(fR11)-1
    if abs(fR11(iG2,1)-fR11(iG2+1,1))<maxpick&&abs(fR11(iG2,2)-fR11(iG2+1,2))<maxpick
        fR11(iG2+1,:)=0;
    end
end
[LG1,WG1]=find(fR11(:,1)==0);
fR11(LG1,:)=[];
pickpoint=fR11;
for ip1=1:length(pickpoint)
    if R(pickpoint(ip1,1),pickpoint(ip1,2))>=150
        pickpoint(ip1,:)=0;
    elseif G(pickpoint(ip1,1),pickpoint(ip1,2))>=200
         pickpoint(ip1,:)=0;
     elseif B(pickpoint(ip1,1),pickpoint(ip1,2))>=220
          pickpoint(ip1,:)=0;
    elseif pickpoint(ip1,:)<0.4/dt
        pickpoint(ip1,:)=0;
    end
end
[Lpickpoint,Wpickpoint]=find(pickpoint(:,1)==0);
pickpoint(Lpickpoint,:)=[];
fpickpoint=pickpoint(:,1:2);
fp=1;
 tfpickpoint=zeros(1,2);
[Lfpickpoint,Wpickpoint]=size(fpickpoint);
for ifpickpoint=1:Lfpickpoint
    fpickpoint(ifpickpoint,3)=vmin+dv*fpickpoint(ifpickpoint,2);
      fpickpoint(ifpickpoint,4)=0.02*fpickpoint(ifpickpoint,3);
end
for ifpickpoint=1:Lfpickpoint-1
    if abs(fpickpoint(ifpickpoint,3)-fpickpoint(ifpickpoint+1,3))<fpickpoint(ifpickpoint,4)
        if abs(fpickpoint(ifpickpoint,1)-ms(fpickpoint(ifpickpoint+1,1)))<0.15/dt;
        if ms(fpickpoint(ifpickpoint,1),fpickpoint(ifpickpoint,2))<ms(fpickpoint(ifpickpoint+1,1),fpickpoint(ifpickpoint+1,2))
            tfpickpoint(fp,:)=fpickpoint(ifpickpoint,1:2);
            fp=fp+1;
        else
            tfpickpoint(fp,:)=fpickpoint(ifpickpoint+1,1:2);
            fp=fp+1;
        end
    end
end
[Ltfpickpoint,Wtfpickpoint]=size(tfpickpoint);
for itfpickpoint=1:Ltfpickpoint
    if tfpickpoint(itfpickpoint,1)==0
        continue
    else
    R(tfpickpoint(itfpickpoint,1),tfpickpoint(itfpickpoint,2))=255;
    G(tfpickpoint(itfpickpoint,1),tfpickpoint(itfpickpoint,2))=255;
    B(tfpickpoint(itfpickpoint,1),tfpickpoint(itfpickpoint,2))=255;
    end
end
for fifR11=1:length(fR11)
    for itfpickpoint=1:Ltfpickpoint
        if fR11(fifR11,1)==tfpickpoint(itfpickpoint,1)&&fR11(fifR11,2)==tfpickpoint(itfpickpoint,2);
            fR11(fifR11,:)=0;
        end
    end
end
[LG1,WG1]=find(fR11(:,1)==0);
fR11(LG1,:)=[];
for iffr11=1:length(fR11)
fR11(iffr11,4)=R(fR11(iffr11,1),fR11(iffr11,2));
fR11(iffr11,5)=G(fR11(iffr11,1),fR11(iffr11,2));
fR11(iffr11,6)=B(fR11(iffr11,1),fR11(iffr11,2));
end
RGB=cat(3,R,G,B); 
for ifR11=1:length(fR11)
    gg(ifR11,:)=fR11(ifR11,:);
end

end
                          end

%% Automatic velocity analysis based on RGB space mapping 
% algorithm
%%Code by:  Junming Zhang, Bin Hu, Deli Wang and Xiangbo Gong
%%email: zhangjm19@mails.jlu.edu.cn  
%% input data
clear ;
  load('mymap.mat');%colormap
load('data.mat');% full wavefield data in shot gather
load('absorb.mat')% absorb boundary data
%% multiple prediction
data2=data;
xmute1=[1,22,43,64,87,106,125];
tmute1=[186,200,237,287,351,405,463];
ntaper=1;
mutype='SRM'; %% free-surface multiple prediction
[SRME] = MultiplePrediction(data2,xmute1,tmute1,ntaper,mutype);
 xmute2=[1,26,54,75,99,125];
 tmute2=[646,649,661,681,708,748];
data3=Mute(data2,xmute2,tmute2,ntaper,1);
mutype='IM'; %% internal multiple prediction 
[IM1]=MultiplePrediction(data3,xmute1,tmute1,ntaper,mutype);
data4= Mute(data3,xmute1,tmute1,ntaper,0);
xmute3=[1,12,35,55,82,103,125];
tmute3=[309,313,324,344,388,430,486];
[IM2]=MultiplePrediction(data4,xmute3,tmute3,ntaper,mutype);
IM=IM1+IM2;
FM=SRME+IM;
Suminput=mean(mean(abs(data2),3),2);
Sumoutput=mean(mean(abs(FM),3),2);
[pks,locs]=findpeaks(Sumoutput,'minpeakdistance',10);
aa=pks(1);
bb=locs(1);
scale=mean(Suminput(bb:end))/mean(Sumoutput(bb:end));
FM=FM*scale;
[FM1,FM2,FM3]=size(FM);
%%Convert dataset to cmp set
for idata=1:FM3
data1(:,:,idata)=data(:,:,1);
absorb1(:,:,idata)=absorb(:,:,1);
end
for iFM3=1:FM3
    datafull(:,FM2*iFM3-FM2+1:FM2*iFM3)=data1(:,:,iFM3);
    FMbw(:,FM2*iFM3-FM2+1:FM2*iFM3)=FM(:,:,iFM3);
    absorbfull(:,FM2*iFM3-FM2+1:FM2*iFM3)=absorb1(:,:,iFM3);
end
[datacmp]=window(datafull,FM2,20,20,FM3,1);
 dataslice=datacmp(:,:,1);
 Fullmutipleaw=window(FMbw,FM2,20,20,FM3,1);
Fullmutiple=Fullmutipleaw(:,:,1);
[absorbcmp]=window(absorbfull,FM2,20,20,FM3,1);
 absorbslice=absorbcmp(:,:,1);
%%Velocity Spectrum
offset=-2480:40:2480;%offset
offsetkm=offset/1000;
dt=0.004;%time sample
N1=8;
 vmin=1000;%Minimum scan velocity
 vmax=2500;%maxium scan velocity
 nv=100;%scan number
 dv=(vmax-vmin)/(nv-1);
 vel=vmin:dv:vmax;
 maxpick=6;
  s1=velan(dataslice,dt,offset,vmin,vmax,nv,1,3,[3 21]);
   s2=velan(Fullmutiple,dt,offset,vmin,vmax,nv,1,3,[3 21]);
   c=velan(absorbslice,dt,offset,vmin,vmax,nv,1,3,[3 21]);
%    sori=velan(bb,dt,offset,vmin,vmax,nv,1,3,[3 21]);

%%RGB mapping
[nt,nh]=size(s1);
 ms=medfilt3(s1,[12,5]);
  mp=medfilt3(s2,[12,5]);
[simi]=Localsimilarity(s1,s2,8,2,0.01,0.001,10,ls);% calculate local similarity
border='Y'; %% Y,N Is it necessary to exclude energy groups with  abnormally high and low velocities
limit='Y'; % Y,N,Whether the selection limit of each feature is needed to constrain the primary pickup
pl=120;%limit for multiple predictability
pv=160;%limit for velocity
pa=180;%limit for amplitude
[shape,point] = RGBmapping(ms,simi,N1,dt,vmin,vmax,nv,maxpick,border);
[RGB,velocityline1,pickpoint]=Pick(ms,shape,point,dt,vmin,dv,limit,pl,pv,pa);
                                  
% %manual pick
c=medfilt3(c,[12,5]);
vstart=(1400-vmin)/dv;
x1=[1,166,293,392,626,767,887,992,1092,1250];
    y1=[vstart,35,50,65,56,57,59,62,65,78];
    manpickpoint=[x1;y1];
    manpickpoint=manpickpoint';
    manpickpoint(:,1)=manpickpoint(:,1).*dt;
    manpickpoint(:,2)=manpickpoint(:,2).*dv+vmin;
 sl2=interp1(x1,y1,1:nt,'cubic');
velocityline2=triangle(20,1,nt,sl2);
% %variational method
cx1=[1,nt];
cy1=[0.2*nh,0.62*nh];
cx2=[1,nt];
cy2=[0.65*nh,0.9*nh];
    cut1=polyfit(cx1,cy1,1);
    cut2=polyfit(cx2,cy2,1);
    line1=interp1(cx1,cy1,1:nt);
    line2=interp1(cx1,cy2,1:nt);
 mx=triangle2(2,10,nh,nt,s1);
d=mx;
  for id1=1:nt 
         d(id1,1:floor(id1*cut1(:,1)+cut1(:,2)))=0;
         d(id1,floor(id1*cut2(:,1)+cut2(:,2)):nh)=0;
  end       
   twindow=20;
   velpick=fomel(d,twindow,vel,dv);

velpick=triangle(15,1,nt,velpick);
pv=vmin:dv:vmax;
pt=dt:dt:nt*dt;
velocityline1=velocityline1.*dv+vmin;
velocityline2=velocityline2.*dv+vmin;
velocityline3=velpick.*dv+vmin;
line1=line1.*dv+vmin;
line2=line2.*dv+vmin;
fosi=15;
[La,Wa]=size(dataslice);
for ia=1:Wa
    fidataslice(:,ia)=dataslice(:,ia)/norm(dataslice(:,ia));
    fiabsorbslice(:,ia)=absorbslice(:,ia)/norm(absorbslice(:,ia));
end
figure(1)
imagesc(offsetkm,pt,clip(fidataslice,0.1))
set(gcf,'color','white','unit','normalized','position',[0.2,-0.3,0.28,1.2])
colormap gray;
title('full wavefield data','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi);
xlabel('offset / km','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi)
ylabel('Time / s','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi)
figure(2)
imagesc(offsetkm,pt,clip(fiabsorbslice,0.1))
set(gcf,'color','white','unit','normalized','position',[0.2,-0.3,0.28,1.2])
colormap gray;
title('predicted multiple','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi);
xlabel('offset / km','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi)
ylabel('Time / s','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi)

figure(3)
imagesc(pv,pt,clip(ms,0.03))
hold on
set(gcf,'color','white','unit','normalized','position',[0.2,-0.3,0.28,1.2])
colormap(newColorMap);
title('full wavefield data','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi);
xlabel('Velocity / m/s','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi)
ylabel('Time / s','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi)

figure(4)
imagesc(pv,pt,clip(mp,0.03))
hold on
set(gcf,'color','white','unit','normalized','position',[0.2,-0.3,0.28,1.2])
colormap(newColorMap);
title('predicted multiple','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi);
xlabel('Velocity / m/s','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi)
ylabel('Time / s','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi)

figure(3)
imagesc(pv,pt,clip(c,0.01))
hold on
plot(velocityline2,pt,'color',[1,0,0],'linewidth',1.5);
set(gcf,'color','white','unit','normalized','position',[0.2,-0.3,0.28,1.2])
colormap(newColorMap);
title('absorbing boundary velocity','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi);
xlabel('Velocity / m/s','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi)
ylabel('Time / s','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi)

figure(4)
imagesc(pv,pt,clip(ms,0.03))
hold on
plot(velocityline3,pt,'color',[1,0,0],'linewidth',1.5);
plot(line1,pt,'color',[1,1,1],'linewidth',1.5);
plot(line2,pt,'color',[1,1,1],'linewidth',1.5);
set(gcf,'color','white','unit','normalized','position',[0.2,-0.3,0.28,1.2])
colormap(newColorMap);
title('variational velocity','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi);
xlabel('Velocity / m/s','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi)
ylabel('Time / s','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi)

figure(5)
plot(pt,velocityline1,'--','color',[1,0,0],'linewidth',2);
axis([0 inf,vmin vmax])
 hold on
plot(manpickpoint(:,1),manpickpoint(:,2),'-.','color',[0,0,1],'linewidth',4,'MarkerSize',10);
plot(pt,velocityline3,'color',[0,1,0],'linewidth',1.5);
 set(gcf,'color','white','unit','normalized','position',[0.2,-0.3,1.1,0.65])
ylabel('Velocity / m/s','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi)
xlabel('Time / s','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi)
legend('RGB line','manpick','tridational');
l1=legend('RGB space','manualpick','variational integral');
set(l1,'Fontname', 'Times New Roman','FontWeight','bold','FontSize',10)

figure(6)
imagesc(pv,pt,RGB)
hold on
plot(velocityline1,pt,'color',[1,0,0],'linewidth',1.5);
set(gcf,'color','white','unit','normalized','position',[0.2,-0.3,0.28,1.2])
colormap(newColorMap);
title('RGB shapping velocity','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi);
xlabel('Velocity / m/s','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi)
ylabel('Time / s','FontName','Times New Roman','FontWeight','Bold','FontSize',fosi)

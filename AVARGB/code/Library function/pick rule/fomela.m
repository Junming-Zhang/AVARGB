function [velpick,vline] = fomela(s1,twindow,vel,dv)
%  twindow=70;
% %    velpick=fomela(ms,twindow,vel,dv);
%    s1=mx;
%    
% s1=pure;
thr1=mean(mean(s1));
[L1]=find(s1>thr1);
thr2=mean(mean(s1(L1)));
aa=find(abs(s1)<thr2);
 s1(aa)=0;

vmin=vel(:,1);
sbackup=s1;
[nt,nh]=size(sbackup);
% twindow=40;
% vel=1000:20:2980;
% lamda=10000;
sbackup=[zeros(twindow,100);sbackup;zeros(twindow,100)];
 ij=1;
for ii=1:nt
    iii=ii+twindow;
    data=sbackup(iii-twindow:iii+twindow,:);
    Vel=repmat(vel,2*twindow+1,1);
    [aa SUM]=max(sum(data,1));
    if SUM<(1430-vmin)/dv
        vinter(ij,1)=ii;
        ij=ij+1;
        SUM=0;
    end
     vline(ii,1)=ii;
     vline(ii,2)=SUM;
    end
%    for i=1:length(vinter)
    
      vline(vinter(:,1),:)=[];  
      vline(1,:)=[1,430/dv];
      vline=[vline;[1250,77]];
      
%        end
%    end
x=vline(:,1);
y=vline(:,2);
    sl2=interp1(x,y,1:1250,'pchip');
 
    velpick=sl2;
 end
 function [velpick] = fomel(s1,twindow,vel,dv)

vmin=vel(:,1);
sbackup=s1;
[nt,nh]=size(sbackup);
% twindow=40;
% vel=1000:20:2980;
% lamda=10000;
sbackup=[zeros(twindow,100);sbackup;zeros(twindow,100)];
for ii=1:nt
    iii=ii+twindow;
    data=sbackup(iii-twindow:iii+twindow,:);
    Vel=repmat(vel,2*twindow+1,1);
    [aa SUM]=max(sum(data,1));
    if SUM<(1430-vmin)/dv
        SUM=(1430-vmin)/dv;
    end
    if max(sum(data,1))==0&&ii>1
        SUM=velpick(ii-1);
    end
    velpick(ii)=SUM;
end

% velpick1=[ones(1,twindow)*velpick(1),velpick,ones(1,twindow)*velpick(1250)];
% 
% for ii=1:1250
%     iii=ii+twindow;
%     data=sbackup(iii-twindow:iii+twindow,:);
%     vel=velpick1(1,iii-twindow:iii+twindow)';
%     Vel=repmat(vel,1,100);
%     [aa SUM]=max(sum(exp(data).*(Vel.^2+lamda.^2).^0.5,1));
%     velpick11(ii)=SUM;
% end

% figure,imagesc(s1)
% hold on 
% plot(velpick,1:1250)
% plot(velpick11,1:1250)
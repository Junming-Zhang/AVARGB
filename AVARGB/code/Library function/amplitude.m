function [amplitudemax,amplitudemin] = amplitude(s1,twindow)

% vmin=vel(:,1);
sbackup=s1;
[nt,nh]=size(sbackup);
% twindow=40;
% vel=1000:20:2980;
% lamda=10000;
sbackup=[zeros(twindow,100);sbackup;zeros(twindow,100)];
for ii=1:nt
    iii=ii+twindow;
    data=sbackup(iii-twindow:iii+twindow,:);
    [Ldata,Wdata]=find(data>0);
    Lpick=[Ldata,Wdata];
    [Lpick1,Wpick1]=size(Lpick);
    for iLpick=1:Lpick1
    datapick(iLpick)=data(Lpick(iLpick,1),Lpick(iLpick,2));
    end
   amax=max(datapick);
   amin=min(datapick);
    amplitudemax(ii)=amax;
    amplitudemin(ii)=amin;
end
end

function [dataout,mask] = Mute(datain,xmute,tmute,ntaper,mode)
% mode=0,1 up or down
len=size(xmute,2);
 [n1 n2]=size(datain);
 xmute(1)=min([1,xmute(1)]);
 xmute(len)=max([xmute(len),n2]);
 xm=1:n2;
 tm=floor(interp1(xmute,tmute,xm,'linear'));   
if (mode==1)
    zone=(ntaper:-1:1)*(1/ntaper);
    data=ones(n1+ntaper,n2);
for ii=1:n2
   data(tm(ii)+ntaper+1:end,ii)=0;
   data(tm(ii)+1:tm(ii)+ntaper,ii)=zone;
end
    mask=data(ntaper+1:end,:);  
    dataout=datain.*mask;
else   
   zone=(1:ntaper)*(1/ntaper);
   data=ones(n1+ntaper,n2);
for ii=1:n2
   data(1:tm(ii)-ntaper,ii)=0;
   data(tm(ii)-ntaper+1:tm(ii),ii)=zone;
end
    mask=data(ntaper+1:end,:);
    dataout=datain.*mask;
end
end


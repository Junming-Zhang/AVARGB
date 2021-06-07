   function [data,CDP] = window(a,s,dp,ds,nh,dcdp)
  %dp=distance shot;
  %ds=distance point
  %s=number of shot
 d=2*(ds/dp);
for ii=1:s
    for jj=1:nh
%         recindex((ii-1)*12+jj)=50+(jj-1)*5+(ii-1)*10;
%         shotindex((ii-1)*12+jj)=50+(ii-1)*10;
         cdpindex((ii-1)*nh+jj)=jj+d*(ii-1);
    end
end
mcdpindex=max(cdpindex);
for ii=1:mcdpindex
aa=find(cdpindex==ii);
disp(['cdp',num2str(ii),':fold ',num2str(size(aa,2))]);
fcdp(ii,1)=ii;
fcdp(ii,2)=size(aa,2);
end
maxfold=max(fcdp(:,2));
[fullcdp,ff]=find(fcdp(:,2)==maxfold);
% cdpmin=min(fullcdp);
% cdpmax=max(fullcdp);

CDP=fullcdp(1:dcdp:end,:);
[lCDP,WCDP]=size(CDP);


for ig=1:lCDP%1:13
%     disp(['ig=',ig]);
    aa=find(cdpindex==CDP(ig));
    
     data2=a(:,aa);
% data=a(:,aa);
    data3=fliplr(data2);
    data4=data3(:,2:length(aa));
     data5=[data2,data4];
       data(:,:,ig)=data5;
end

%     data(:,:,ii)=[a(:,aa) fliplr(a(:,aa(1:length(aa-1))))];
%     VelocityAnalysis(data, t, x, vmin, vmax, vstep);
end

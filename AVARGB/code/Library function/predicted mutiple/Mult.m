clear 
load data2
[nt ntrace]=size(data);
data=data(:,1:2:end);
xmute=[1,52,127,244,272,290,301];
tmute=[334,342,387,499,530,556,561];
ntaper=10;

mutype='SRM'; %% SRM,IM,FM 
[Srme] = MultiplePrediction(data,xmute,tmute,ntaper,mutype);
mutype='IM'; %% SRM,IM,FM 
[Cfp] = MultiplePrediction(data,xmute,tmute,ntaper,mutype);
mutype='FM'; %% SRM,IM,FM 
[Fm] = MultiplePrediction(data,xmute,tmute,ntaper,mutype);


figure(2),
subplot(141),imagesc(data);colormap gray
subplot(142),imagesc(Srme(:,:,1));colormap gray
subplot(143),imagesc(Cfp(:,:,1));colormap gray
subplot(144),imagesc(Fm(:,:,1));colormap gray



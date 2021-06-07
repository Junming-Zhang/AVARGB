function [output] = MultiplePrediction(input,xmute,tmute,ntaper,mutype)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function used for obtaining the preditced multiple
% Code by:  Junming Zhang, Bin Hu, Deli Wang and Xiangbo Gong 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xmute: The index of muting zone
% tmute: The index of muting zone
% ntaper: The length of taper operater
% mutype: The type of predicted multiple 
% SRM: Predict surfaced related multiple only
% IM: Predict internal multiple only
% FM: Predict fullwavefield multiple 

[nt ntrace]=size(input);
srme=zeros(nt,ntrace,ntrace);
cfp=zeros(nt,ntrace,ntrace);
%   Surface Related Multiple Prediction


if (strcmp(mutype,'SRM'))
[data] = DataRegularization(input);
output = SRME(data,data);
[srme] = DataRegularization_back(output);
clear output data
end

%   Surface Related Multiple Prediction
if (strcmp(mutype,'IM'))
[datadown,mask] = Mute(input,xmute,tmute,ntaper,0);
dataup=input-datadown;
[dataup] = DataRegularization(dataup);
[datadown] = DataRegularization(datadown);
dataup = SRME(datadown,dataup,'Corr');
dataup = SRME(dataup,datadown);
[cfp] = DataRegularization_back(dataup);
clear dataup datadown
end
output=srme+cfp;

end


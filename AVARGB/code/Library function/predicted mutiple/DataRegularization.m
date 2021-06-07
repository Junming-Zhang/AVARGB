function [output] = DataRegularization(input)
%DataRegularization
[nt ntrace]=size(input);
input1=fliplr(input);
data=[input1(:,1:ntrace-1),input];
ns=ntrace;
output=zeros(nt,ntrace,ns);
for ii=1:ns
    ishot=ntrace+1-ii;
    output(:,:,ii)=data(:,ishot:ishot-1+ntrace);
end

end


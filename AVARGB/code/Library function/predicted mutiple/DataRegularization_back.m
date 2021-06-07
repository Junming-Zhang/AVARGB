function [output] = DataRegularization_back(input)
%DataRegularization
[nt ntrace ns]=size(input);

output=zeros(nt,ntrace,ns);
for ii=1:ns
    output(:,:,ii)=input(:,:,1);
end

end


function [ output_output ] = SRME( input_A,input_B, type )
%srme 
%   out=SRME(A,B);
if nargin <3
    type='Conv';
end
n1  = size(input_A,1);
n2  = size(input_A,2);
n3 = size(input_A,3);
nt_conv = 2*n1; % time sample length of padded kernel for convolution 
Suminput=mean(mean(abs(input_A),3),2);
INput=fft([zeros(n1,n2,n3);input_A],[],1);
clear input_A;
INput = permute(INput, [2 3 1]); 
INput2=fft([zeros(n1,n2,n3);input_B],[],1);
clear input_B;
INput2 = permute(INput2, [2 3 1]);  
if strcmp(type,'Corr')
    for kkk=1:nt_conv
        INput(:,:,kkk) = -1.* (INput(:,:,kkk)*(INput2(:,:,kkk))');
    end
 clear INput2;
INput = permute(INput,[3 1 2]);
INput= real(ifft(INput,[],1));
output_output=INput(1:n1,:,:);
else
    for kkk=1:nt_conv
        INput(:,:,kkk) = -1.* (INput(:,:,kkk)*(INput2(:,:,kkk)).');
    end
 clear INput2;
INput = permute(INput,[3 1 2]);
INput= real(ifft(INput,[],1));
output_output=INput(1:n1,:,:);
end
Sumoutput=mean(mean(abs(output_output),3),2);
[pks,locs]=findpeaks(Sumoutput,'minpeakdistance',10);
aa=pks(1);
bb=locs(1);
scale=mean(Suminput(bb:end))/mean(Sumoutput(bb:end));
output_output=output_output*scale;
clear input_A
clear INput
end


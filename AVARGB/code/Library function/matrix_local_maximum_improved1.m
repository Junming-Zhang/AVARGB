function ALocalMax=matrix_local_maximum_improved1(A,C)
% prompt='Input row or column number of local region : ';
numofrl=num2str(C);% char
%%construct matrix
B=nan(size(A,1)+2*str2double(numofrl),size(A,2)+2*str2double(numofrl));
B(str2double(numofrl)+1:end-str2double(numofrl),str2double(numofrl)+1:end-str2double(numofrl))=A;
ALocalMaxL=zeros(size(B));

%%find local maximum
for ii=str2double(numofrl)+1:size(A,1)+str2double(numofrl)
    for jj=str2double(numofrl)+1:size(A,2)+str2double(numofrl)
       if B(ii,jj)==max(max(B(ii-str2double(numofrl):ii+str2double(numofrl),jj-str2double(numofrl):jj+str2double(numofrl))))
           ALocalMaxL(ii,jj)=max(max(B(ii-str2double(numofrl):ii+str2double(numofrl),jj-str2double(numofrl):jj+str2double(numofrl))));
       end
    end
end


ALocalMax=ALocalMaxL(str2double(numofrl)+1:end-str2double(numofrl),str2double(numofrl)+1:end-str2double(numofrl));
end
 function [w] = AHP(A)
 [L,W]=size(A);
for i=1:W
    weight(:,i)=A./A(:,i);
end
[x,y]=eig(weight);
lamda=max(diag(y)); 
num=find(diag(y)==lamda); 
w0=x(:,num)/sum(x(:,num)); 
ri=[0,0,0.58,0.90,1.12,1.24,1.32,1.41,1.45]; %一致性指标 
cr0=[(lamda-3)/(3-1)]/ri(3) ;
if abs(cr0)<0.1
x1=sum(x(:,1));
w=x(:,1)./x1;
w=w';
else 
 disp('Error:consistency test failed:');
end


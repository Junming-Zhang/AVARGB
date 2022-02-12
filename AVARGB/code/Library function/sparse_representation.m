function [ Xout ] = sparse_representation(X1,scales,dipc,per)
%per影响滤波效果，越小效果越好，但需要注意边界位置
%scales似乎影响不大，但需要进一步验证
%
  sl2dp=SLgetShearletSystem2D(0,size(X1,1),size(X1,2),scales,dipc);
  coeffsp=SLsheardec2D(X1,sl2dp);
  aa=sort(abs(coeffsp(:)),'descend');
  nlen=size(aa,1);
  thresh=aa(floor(nlen*per));
  aa=find(abs(coeffsp)<thresh);
   coeffsp(aa)=0;
  Xout = SLshearrec2D(coeffsp,sl2dp);
end


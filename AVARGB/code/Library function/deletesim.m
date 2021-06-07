function [yq11] = deletesim(yq11)
yq11=sortrows(yq11,1);
ua1=unique(yq11(:,1));   %找出第一列的惟一值
[h s]=hist(yq11(:,1),ua1); %以这个惟一值做直方图,h是数量，s是相应取值
sh1=s(h~=1);           %找出s中不是1的，就是A(:,1)中重复的
ind=ismember(yq11(:,1),sh1); %A(:,1)中包含sh1的位置
 yq11(ind,:)=[];       %这个位置的这一行删除掉
end
function [yq11] = deletesim(yq11)
yq11=sortrows(yq11,1);
ua1=unique(yq11(:,1));   %�ҳ���һ�е�Ωһֵ
[h s]=hist(yq11(:,1),ua1); %�����Ωһֵ��ֱ��ͼ,h��������s����Ӧȡֵ
sh1=s(h~=1);           %�ҳ�s�в���1�ģ�����A(:,1)���ظ���
ind=ismember(yq11(:,1),sh1); %A(:,1)�а���sh1��λ��
 yq11(ind,:)=[];       %���λ�õ���һ��ɾ����
end
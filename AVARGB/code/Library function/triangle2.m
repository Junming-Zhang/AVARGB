function [ vv ] = triangle2(rect1, rect2, n1, n2, uu)
%vv=triangle2(3,3,6,4,a,b);
 for m=1:n1;
          d(:,m) =  triangle( rect2, n1, n2, uu(:,m));%�ݷ����˲�
 end;

         %call triangle( rect2, n1, n2, uu(i1,1), ss(i1,1))

          ;
 for m=1:n2;
 vv(m,:) =   triangle( rect1,  n2, n1,d(m,:));
 end;

           %  call triangle( rect1,  1, n1, ss(1,i2), vv(1,i2))
%  subroutine triangle2( rect1, rect2, n1, n2, uu, vv)
% integer i1,i2,        rect1�˲������, rect2�˲�������, n1���ݳ���, n2���ݿ��
% real uu(n1,n2), vv(n1,n2)
% temporary real ss(n1,n2)
% do i1= 1, n1
%                 call triangle( rect2, n1, n2, uu(i1,1), ss(i1,1)һ���������)
% do i2= 1, n2
%                 call triangle( rect1,  1, n1, ss(1,i2), vv(1,i2)�����������)
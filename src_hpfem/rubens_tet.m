xyztet=[-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1];
x=xyztet(:,1);
y=xyztet(:,2);
z=xyztet(:,3);


V01=1/6*(  x(2)*(y(3)*z(4)-y(4)*z(3))+ x(3)*(y(4)*z(2)-y(2)*z(4))+ x(4)*(y(2)*z(3)-y(3)*z(2)));

V02=1/6*(  x(1)*(y(4)*z(3)-y(3)*z(4))+ x(3)*(y(1)*z(4)-y(4)*z(1))+ x(4)*(y(3)*z(1)-y(1)*z(3)));

V03=1/6*(  x(1)*(y(2)*z(4)-y(4)*z(2))+ x(2)*(y(4)*z(1)-y(1)*z(4))+ x(4)*(y(1)*z(2)-y(2)*z(1)));

V04=1/6*(  x(1)*(y(3)*z(2)-y(2)*z(3))+ x(2)*(y(1)*z(3)-y(3)*z(1))+ x(3)*(y(2)*z(1)-y(1)*z(2)));

vol=V01+V02+V03+V04

xx=zeros(4,4);
yy=zeros(4,4);
zz=zeros(4,4);

for i=1:4
for j=1:4
xx(i,j)=x(i)-x(j);
yy(i,j)=y(i)-y(j);
zz(i,j)=z(i)-z(j);
end
end


mat=1/6/vol*[ 6*V01 yy(4,2)*zz(3,2)-yy(3,2)*zz(4,2)  xx(3,2)*zz(4,2)-xx(4,2)*zz(3,2)  xx(4,2)*yy(3,2)-xx(3,2)*yy(4,2);
              6*V02 yy(3,1)*zz(4,3)-yy(3,4)*zz(1,3)  xx(4,3)*zz(3,1)-xx(1,3)*zz(3,4)  xx(3,1)*yy(4,3)-xx(3,4)*yy(1,3);
	          6*V03 yy(2,4)*zz(1,4)-yy(1,4)*zz(2,4)  xx(1,4)*zz(2,4)-xx(2,4)*zz(1,4)  xx(2,4)*yy(1,4)-xx(1,4)*yy(2,4);
	          6*V04 yy(1,3)*zz(2,1)-yy(1,2)*zz(3,1)  xx(2,1)*zz(1,3)-xx(3,1)*zz(1,2)  xx(1,3)*yy(2,1)-xx(1,2)*yy(3,1)]
	      
	      
% Now Check
for i=1:4
    vec=[1 xyztet(i,:)]';
    mat*vec
    pause
end

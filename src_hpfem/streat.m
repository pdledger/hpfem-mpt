function [unkz,unksid,unkfatp1,unkfatp2,unkfatp3,unkint,sol]=streat(known,npec,sol,...
    nunkt,unksid,order,nside,nface,nelem,nrhs,unkfatp1,unkfatp2,unkfatp3,unkz,unkint)


for i=1:nside
    if unkz(i)<0
        unkz(i)=abs(unkz(i))+nunkt;
    end
end

for i=1:nside
    for j=1:order
        if unksid(i,j)<0
            unksid(i,j)=abs(unksid(i,j))+nunkt;
        end
    end
end

for i=1:nface
    for j=1:(order*order-order)/2
        if unkfatp1(i,j)<0
            unkfatp1(i,j)=abs(unkfatp1(i,j))+nunkt;
        end
    end
end

for i=1:nface
    for j=1:(order*order-order)/2
        if(unkfatp2(i,j)<0)
            unkfatp2(i,j)=abs(unkfatp2(i,j))+nunkt;
        end
    end
end

for i=1:nface
    for j=1:order-1
        if unkfatp3(i,j)<0
            unkfatp3(i,j)=abs(unkfatp3(i,j))+nunkt;
        end
    end
end

sol(1+nunkt:npec+nunkt,1:nrhs)=known(1:npec,1:nrhs);

nunktotal=size(sol,1);
% correct for zero terms fix all zero values to be at nunktotal+1
sol(nunktotal+1,1:nrhs)=zeros(1,nrhs);

for i=1:nside
    if unkz(i)==0
        unkz(i)=nunktotal+1;
    end
end

for i=1:nside
    for j=1:order
        if unksid(i,j)==0
            unksid(i,j)=nunktotal+1;
        end
    end
end

for i=1:nface
    for j=1:(order*order-order)/2
        if unkfatp1(i,j)==0
            unkfatp1(i,j)=nunktotal+1;
        end
    end
end

for i=1:nface
    for j=1:(order*order-order)/2
        if(unkfatp2(i,j)==0)
            unkfatp2(i,j)=nunktotal+1;
        end
    end
end

for i=1:nface
    for j=1:order-1
        if unkfatp3(i,j)==0
            unkfatp3(i,j)=nunktotal+1;
        end
    end
end

if order>=3
for i=1:nelem
    for j=1:(order-2)*(order-1)*(order+1)/2
        if unkint(i,j)==0
            unkint(i,j)=nunktotal+1;
        end
    end
end
end
% Function for reading a basic ANSYS mesh

function [quadlin,intma,nelem,coord,mat,npoin,cintma]= ansysmeshinfo();
% read in from an-Ansys mesh file.
nboun=[];
A=load('connectivity-1.txt');
B=load('nodes-1.txt');
disp('done')
intma=A(:,2:5);
[nelem dum]=size(intma);
cintma=A(:,6:11);
mat=A(:,12);
coord=B(:,2:4);
[npoin dum]=size(coord);
bcosurf=[];
nsurf=[];
quadlin=2;

dintma=cintma;
cintma(:,6-4)=dintma(:,7-4);
cintma(:,7-4)=dintma(:,8-4);
cintma(:,8-4)=dintma(:,6-4);

tol=1e-10;
ndup=0;
flag=zeros(npoin,1);
for i=1:npoin
    for j=[i+1:1:npoin] % only check subsequent occurances of this node.
        if norm(coord(i,:)-coord(j,:))<=tol & flag(j)==0
            % this is a duplicate node
            flag(j)=i;
            ndup=ndup+1;
        end
    end
end


disp(['We have found',num2str(ndup),'nodes'])
% remove duplicates
for i=1:nelem
    for j=1:4
        if flag(intma(i,j))~=0
            intma(i,j)=flag(intma(i,j));
        end
    end
    for j=1:6
        if flag(cintma(i,j))~=0
            cintma(i,j)=flag(cintma(i,j));
        end
    end
end
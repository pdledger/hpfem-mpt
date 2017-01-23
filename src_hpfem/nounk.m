% Build the unknown numbering on edges, faces and element interiors
% assuming a uniform order elements throughout the mesh.

% Dirichlet DOFs are numbered to be negative.

% Edge, Face and Interior numbering are seperated in to the blocks
% Zero block
% High order edge gradient block
% High order face gradient block
% High order face non-gradient block
% High order interior gradient block
% High order interior non-gradient block


function [unkz,unksid,unkfatp1,unkfatp2,unkfatp3,unkint,nef,npec,nunk,nunkt,...
    ordersid,orderfa,help1]=nounk(nelem,nside,nface,cond,glob,globfa,order,...
    orderel,bcedge,matc)

% unksid = unknown numbers on element edges
% unkint = unknown numbers on element interiors
% unkfac = unknown numbers on element faces

% set up help array
help1 = zeros(1,nface);
for i=1:nelem
    for j=1:4
        if cond(i,j)==1 || cond(i,j)==2 || cond(i,j)==4
            help1(globfa(i,j))=cond(i,j);
        end
    end
end

% use the minimum rule to determine p on edges and faces
ordersid = 10000*ones(1,nside);
orderfa = 10000*ones(1,nface);
for i = 1:nelem
    for j=1:6
        if ordersid(glob(i,j))>orderel(i)
            ordersid(glob(i,j)) = orderel(i);
        end
    end
    for j=1:4
        if orderfa(globfa(i,j))>orderel(i)
            orderfa(globfa(i,j)) = orderel(i);
        end
    end
end

% flag edges and faces in Region 1 so as to skip gradients
help2 = zeros(1,nside);
help3 = zeros(1,nface);

% Flag only those edges in the interior of free space
for i=1:nelem
    if matc(i)==0
        for j=1:6
            help2(glob(i,j))=1;
        end
    end
end
for i=1:nelem
    if matc(i)==1
        for j=1:6
            help2(glob(i,j))=0;
        end
    end
end

% Flag only those faces in the interior of free space
for i=1:nelem
    if matc(i)==0
        for j=1:4
            help3(globfa(i,j))=1;
        end
    end
end
for i=1:nelem
    if matc(i)==1
        for j=1:4
            help3(globfa(i,j))=0;
        end
    end
end




unkz = zeros(nside,1);
unksid = zeros(nside,order);

npec=0;
nunk=0;

% treat 1,2,4 Dirichlet (2 posibly non-zero)
% treat 3 Neumann
for i = 1:nside
    if bcedge(i)~=1 && bcedge(i)~=2 && bcedge(i)~=4
        nunk = nunk +1;
        unkz(i) = nunk;
    elseif bcedge(i) ==2
        npec = npec+1;
        unkz(i) = -1*npec;
    end
end

% higher order edges
if order >=1
    for i=1:nside
%        if help2(i)==0
            if bcedge(i)~=1 && bcedge(i)~=2 && bcedge(i)~=4 && help2(i)==0
                for j = 1:ordersid(i)
                    nunk = nunk+1;
                    unksid(i,j) = nunk;
                end
            elseif bcedge(i)==2
                for j=1:order
                    npec = npec+1;
                    unksid(i,j) = -1*npec;
                end
            end
        end
%    end
end
display(['The number of edge based unknowns is ',num2str(nunk)])

% Face functions
unkfatp1=[];
unkfatp2=[];
unkfatp3=[];
if order >=2
    % zero unkfatp1
    unkfatp1 = zeros(nface,(order*order-order)/2);
    for i =1:nface
 %       if help3(i)==0
            if help1(i) ==0 && help3(i)==0
                k = 0;
                for j1 = 0:order-2
                    for j2=0:order-2
                        if (j1+j2)<=(order-2)
                            k=k+1;
                            if (j1+j2)<=(orderfa(i)-2)
                                nunk = nunk+1;
                                unkfatp1(i,k)=nunk;
                            end
                        end
                    end
                end
            elseif help1(i)==2
                k=0;
                for j1=0:order-2
                    for j2=0:order-2
                        if (j1+j2)<=(order-2)
                            k = k+1;
                            npec=npec+1;
                            unkfatp1(i,k)=-1*npec;
                        end
                    end
                end
                
            end
  %      end
    end
    
    % zero unkfatp2
    unkfatp2 = zeros(nface,(order*order-order)/2);
    for i=1:nface
        if help1(i)==0
            k=0;
            for j1=0:order-2
                for j2=0:order-2
                    if (j1+j2)<=(order-2)
                        k=k+1;
                        if (j1+j2)<=(orderfa(i)-2)
                            nunk = nunk+1;
                            unkfatp2(i,k)=nunk;
                        end
                    end
                end
            end
        elseif help1(i)==2
            k=0;
            for j1=0:order-2
                for j2=0:order-2
                    if (j1+j2)<=(order-2)
                        k=k+1;
                        npec=npec+1;
                        unkfatp2(i,k)=-1*npec;
                    end
                end
            end
        end
    end
    
    % zero unkfatp3
    unkfatp3 = zeros(nface,order-1);
    for i = 1:nface
        if help1(i)==0
            for j=1:orderfa(i)-1
                nunk=nunk+1;
                unkfatp3(i,j)= nunk;
            end
        elseif help1(i)==2
            for j = 1:order-1
                npec = npec+1;
                unkfatp3(i,j)=-1*npec;
            end
        end
    end
end
display(['The total number of edge and face unknowns is',num2str(nunk)]);
display(['The total number of PEC edge and face knowns is',num2str(npec)]);

% save number of H curl unknowns
nef = nunk;
nunkt = nunk;

% Interior functions
unkint = [];
if order >=3
    unkint = zeros(nelem,(order-1)*(order-2)*(order+1)/2);
    for i = 1:nelem
        if matc(i)==1
            % gradients
            k=0;
            for ii=0:order-3
                for jj=0:order-3
                    for kk=0:order-3
                        if ii+jj+kk <= order-3
                            k=k+1;
                            nunkt = nunkt+1;
                            unkint(i,k) = nunkt;
                        end
                    end
                end
            end
        end
    end
    % if there are no conductors the above loop may give k=0
    % so as to fix the numbering do the following
     k=0;
            for ii=0:order-3
                for jj=0:order-3
                    for kk=0:order-3
                        if ii+jj+kk <= order-3
                            k=k+1;
                        end
                    end
                end
            end
    % non -gradients (do not distinguish between)
    nintbas=k;
    for i = 1:nelem
        k=nintbas;
        for ii=nintbas+1:(order-2)*(order-1)*(order+1)/2
            nunkt = nunkt+1;
            k=k+1;
            unkint(i,k) = nunkt;
        end
    end
    
    %     for i = 1:nelem
    %         % non-gradients type 2a
    %         k=nintbas;
    %         for ii=0:order-3
    %             for jj=0:order-3
    %                 for kk=0:order-3
    %                     if ii+jj+kk <= order-3
    %                         k=k+1;
    %                         nunkt = nunkt+1;
    %                         unkint(i,k) = nunkt;
    %                     end
    %                 end
    %             end
    %         end
    %         % non-gradients type 2b
    %         for ii=0:order-3
    %             for jj=0:order-3
    %                 for kk=0:order-3
    %                     if ii+jj+kk <= order-3
    %                         k=k+1;
    %                         nunkt = nunkt+1;
    %                         unkint(i,k) = nunkt;
    %                     end
    %                 end
    %             end
    %         end
    %         % non-gradients type 3
    %         for jj=0:order-3
    %             for kk=0:order-3
    %                 if jj+kk <= order-3
    %                     k=k+1;
    %                     nunkt = nunkt+1;
    %                     unkint(i,k) = nunkt;
    %                 end
    %             end
    %         end
    
    %    end
    if k~=(order-2)*(order-1)*(order+1)/2
        disp('wrong number of interior functions created')
    end
    
end


display(['The total number of edge face and interior unkowns is',num2str(nunkt)]);

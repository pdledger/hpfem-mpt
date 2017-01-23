% Function for performing the main matrix assembly and matrices for
% preconditioning operations.

% The following steps are performed-
% Integration points obtained and basis, curl basis and geometry
% interpolatory functions evaluated at integration points.
% Dirichlet DOFs computed (from L2 minmisation) and stored.
%

function [rhs,stmtx,known,ephx1,ephy1,ephz1,ecphx1,ecphy1,ecphz1,ephx2,...
    ephy2,ephz2,ecphx2,ecphy2,ecphz2,intxi,inteta,intzeta,intw,nip,intfxi,...
    intfet,intfw,nipf,premtx,probdata,gphx1,gphy1,gphz1,gphx2,gphy2,gphz2,phh11,phh12]= eigenass(nelem,intma,coord,order,eltype,glob,unksid,...
    nside,cond,globfa,unkint,esizet,nipd,nunkt,nef,edgecof,facecof,...
    gorder,bcedge,npec,orderel,nrhs,unkfatp1,unkfatp2,unkfatp3,mat,unkz,probdata,matc)


display('set mat parameters')
econt=6*(order+1);
maxnip = nipd^3;

lec = zeros(6,gorder,3);
lfc = zeros(4,(gorder*(gorder-1)/2),3);

if(order>=2)
    econt=econt+(12*(order-1))+(4*(order-1)*(order-2));
end

% obtain integration points......over a tetrahedra
[intxi,inteta,intzeta,intw,nip]=intpoints(nipd,maxnip);
if nip~=nipd^3
    error(message('dimension problem for nipd'));
end

% obtain integration points......over a face
[intfxi,intfet,intfw,nipf]=intpointf(nipd,maxnip);

% obtain ph and curl ph at all int points inside tet
[ephx1,ephy1,ephz1,ephx2,ephy2,ephz2,ecphx1,ecphy1,ecphz1,ecphx2,...
    ecphy2,ecphz2,gphx1,gphy1,gphz1,gphx2,gphy2,gphz2,phh11,phh12]= evaluate(nip,intxi,inteta,intzeta,order,esizet,gorder);

display('completed all integration data')

% determine known Dirchlet values from a L2 projection
[known,ephdfx1,ephdfy1,ephdfz1,ephdfx2,ephdfy2,ephdfz2,phh1f1,phh1f2,gphfx1,...
    gphfy1,gphfz1,gphfx2,gphfy2,gphfz2]= pec(nelem,nside,coord,lec,lfc,glob,globfa,intma,...
    edgecof,facecof,gorder,order,unksid,esizet,eltype,bcedge,nipf,...
    intfxi,intfet,intfw,cond,npec,nrhs,unkfatp1,unkfatp2,unkfatp3,unkz,probdata);


% enter main assembly routine
% intialize arrays before proceeding to parrallel region
kc = zeros(esizet,esizet);
kcpre = zeros(esizet,esizet);
gc = zeros(esizet,nrhs);
mdiopole=zeros(3,1);

nmst =0;
I = zeros(nunkt,1);
J = zeros(nunkt,1);
X = zeros(nunkt,1);
Xpre = zeros(nunkt,1);
rhs = zeros(nunkt,nrhs);
sigma=probdata.matr.sigma;
omega=probdata.matr.omega;
muz=probdata.matr.muz;
epz=probdata.matr.epz;
delta=probdata.matr.delta;
mu=probdata.matr.mu;
epl=probdata.matr.epl;
jsrc=probdata.matr.jsrc;

% Determine the number of gradient interior functions.
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
nintbas=k;

% perfrom assembly
disp('Begin main assembly')
for i=1:nelem
    if matc(i)==1
        %     conductor
        kappalow=1i*omega*sigma(mat(i))*muz*delta^2 -omega^2*epl(mat(i))*muz*epz;
    else
        %     regularised free space
        kappalow=omega*sigma(mat(i))*muz*delta^2 -omega^2*epl(mat(i))*muz*epz;
    end
    %     set kappa for the higher order block
    if probdata.sol.regopt ==1
        if matc(i)~=1
            kappa=0;
        else
            kappa=kappalow;
        end
    else
        kappa=kappalow;
    end
    mue=mu(mat(i));
    
    % pelem= order of the element i.
    % use this only as a flag for static condensation
    pelem=orderel(i);
    
    if mod(i,1000)==0
        disp([num2str(i),' Elements of ',num2str(nelem),' Processed']);
    end
    
    if i == nelem
        display('All elem Processed')
    end
    
    % transfer coordinates to local array
    xy = coord(intma(i,1:4),1:3);
    flag=0 ;
    gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
    mycoord = zeros(gesizet,3);
    %transfer coefficents to locations vertices
    mycoord(1:4,1:3) = xy(1:4,1:3);
    lec=zeros(6,gorder,3);
    lfc=zeros(4,gorder*(gorder-1)/2,3);
    if gorder > 0
        for j=1:6
            for p=1:gorder
                for k=1:3
                    lec(j,p,k)=edgecof(glob(i,j),((p-1)*3)+k);
                    if abs(edgecof(glob(i,j),((p-1)*3)+k))>0.00000001
                        flag=1;
                    end
                    mycoord(4+j+6*(p-1),k)=lec(j,p,k);
                end
            end
        end
        
        for j=1:4
            for p=1:gorder*(gorder-1)/2
                for k=1:3
                    lfc(j,p,k)=facecof(globfa(i,j),((p-1)*3)+k);
                    mycoord(4+6*gorder+(j-1)*gorder*(gorder-1)/2+p,k)= lfc(j,p,k);
                end
            end
        end
    end
    
    bhelp=rowfun(unkz,unksid,unkint,order,glob,globfa,...
        i,unkfatp1,unkfatp2,unkfatp3,esizet);
    
    % transfer jsrc to local array
    ljsrc(1:3)=jsrc(mat(i),1:3);
    
    % calculate elemental stiffness, mass and source matrices
    if eltype(i)==1
        stiff= stiffmatrix(esizet,intw,nip,ecphx1,ecphy1,ecphz1,flag,gorder,mue,gphx1,...
            gphy1,gphz1,mycoord);
        
        mass= massmatrix(esizet,intw,nip,ephx1,ephy1,ephz1,flag,gorder,gphx1,...
            gphy1,gphz1,mycoord);
        
        erhs= elemrhs(xy,esizet,intfxi,intfet,intfw,nipf,order,eltype(i),cond,...
            i,lec,lfc,flag,gorder,mue,nrhs,probdata,ephdfx1,ephdfy1,ephdfz1,mat,...
            mycoord,phh1f1,gphfx1,gphfy1,gphfz1)  ;
        
        if matc(i)==1
            [erhs,mdiopole]= jmatrix(esizet,intw,nip,ephx1,ephy1,ephz1,flag...
                ,gorder,ljsrc,erhs,mdiopole,...
                probdata,gphx1,gphy1,gphz1,kappa,mycoord,phh11);
        end
        
    else
        stiff= stiffmatrix(esizet,intw,nip,ecphx2,ecphy2,ecphz2,flag,gorder,mue,...
            gphx2,gphy2,gphz2,mycoord);
        
        mass= massmatrix(esizet,intw,nip,ephx2,ephy2,ephz2,flag,gorder,gphx2,...
            gphy2,gphz2,mycoord);
        
        erhs= elemrhs(xy,esizet,intfxi,intfet,intfw,nipf,order,eltype(i),cond,...
            i,lec,lfc,flag,gorder,mue,nrhs,probdata,ephdfx2,ephdfy2,ephdfz2,mat,...
            mycoord,phh1f2,gphfx2,gphfy2,gphfz2)  ;
        if matc(i)==1
            [erhs,mdiopole]= jmatrix(esizet,intw,nip,ephx2,ephy2,ephz2,flag...
             ,gorder,ljsrc,erhs,mdiopole,...
             probdata,gphx2,gphy2,gphz2,kappa,mycoord,phh12);
        end
        
    end
    
    
    kc(1:6,1:6)=stiff(1:6,1:6)+kappalow*mass(1:6,1:6);
    kcpre(1:6,1:6)=stiff(1:6,1:6)+kappalow*mass(1:6,1:6);
    gc(1:6,1:nrhs)=erhs(1:6,1:nrhs);
    
    %     higher order block
    for p=1:esizet
        if p<= 6
            kc(p,7:esizet)=stiff(p,7:esizet)+kappa*mass(p,7:esizet);
        else
            kc(p,1:esizet)=stiff(p,1:esizet)+kappa*mass(p,1:esizet);
            kcpre(p,1:esizet)=stiff(p,1:esizet)+abs(kappa)*mass(p,1:esizet);
        end
        gc(p,1:nrhs)=erhs(p,1:nrhs);
    end
    
    if probdata.sol.regopt==1
        % high order edges (gradients)
        kcpre(7:6*(order+1),7:6*(order+1))=stiff(7:6*(order+1),7:6*(order+1))+abs(kappa)*mass(7:6*(order+1),7:6*(order+1));
        
        % high order faces (gradients)
        st=6*(order+1)+1;en=6*(order+1)+4*(order*order-order)/2;
        kcpre(st:en,st:en)=stiff(st:en,st:en)+abs(kappa)*mass(st:en,st:en);
        
        % high order faces (non-gradients)
        st=6*(order+1)+4*(order*order-order)/2+1; en=6*(order+1)+4*(2*(order*order-order)/2+(order-1));
        kcpre(st:en,st:en)=stiff(st:en,st:en)+abs(kappalow)*mass(st:en,st:en);
        
        % high order interiors (gradients)
        st=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+1;
        en=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+nintbas;
        kcpre(st:en,st:en)=stiff(st:en,st:en)+abs(kappa)*mass(st:en,st:en);
        
        % high order interiors (non-gradients)
        st=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+nintbas+1;
        en=esizet;
        kcpre(st:en,st:en)=stiff(st:en,st:en)+abs(kappalow)*mass(st:en,st:en);
    else
        % high order edges (gradients)
        kcpre(7:6*(order+1),7:6*(order+1))=stiff(7:6*(order+1),7:6*(order+1))+abs(kappalow)*mass(7:6*(order+1),7:6*(order+1));
        
        % high order faces (gradients)
        st=6*(order+1)+1;en=6*(order+1)+4*(order*order-order)/2;
        kcpre(st:en,st:en)=stiff(st:en,st:en)+abs(kappalow)*mass(st:en,st:en);
        
        % high order faces (non-gradients)
        st=6*(order+1)+4*(order*order-order)/2+1; en=6*(order+1)+4*(2*(order*order-order)/2+(order-1));
        kcpre(st:en,st:en)=stiff(st:en,st:en)+abs(kappalow)*mass(st:en,st:en);
        
        % high order interiors (gradients)
        st=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+1;
        en=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+nintbas;
        kcpre(st:en,st:en)=stiff(st:en,st:en)+abs(kappalow)*mass(st:en,st:en);
        
        % high order interiors (non-gradients)
        st=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+nintbas+1;
        en=esizet;
        kcpre(st:en,st:en)=stiff(st:en,st:en)+abs(kappalow)*mass(st:en,st:en);
    end
    
    % assemble contributions in to global stiffness matrix and right
    % hand-side vector
    for j = 1:esizet%econt
        row = bhelp(j);
        if row>0
            for k = 1:esizet%econt
                col = bhelp(k);
                if col>0
                    nmst=nmst+1;
                    len = length(X);
                    if nmst>len
                        I(2*len)=0;
                        J(2*len)=0;
                        X(2*len)=0;
                        Xpre(2*len)=0;
                    end
                    I(nmst) = row;
                    J(nmst) = col;
                    X(nmst) = kc(j,k);
                    Xpre(nmst) = kcpre(j,k);
                elseif col<0
                    % move Dirihlet columns to right hand side
                    rhs(row,1:nrhs)=rhs(row,1:nrhs)-...
                        ((kc(j,k))*known(abs(col),1:nrhs));
                end
            end
            % Add right hand side source terms
            rhs(row,1:nrhs)=rhs(row,1:nrhs)+gc(j,1:nrhs);
        end
    end
    
end


stmtx = sparse(I(1:nmst),J(1:nmst),X(1:nmst),nunkt,nunkt);
premtx = sparse(I(1:nmst),J(1:nmst),Xpre(1:nmst),nunkt,nunkt);

%clear I
%clear J
%clear X
%clear Xpre

display('completed all assembly for continous dof')



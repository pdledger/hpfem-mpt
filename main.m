% MATLAB implemenetation of the computation of the magnetic PT

clear all, clc, close all
format long
addpath(genpath('./'))

% Options of different geometries
% problem 1 = conducting magnetic sphere
% problem 2 = conducting magnetic spheroid
% problem 3 = nested sphere conducting and magnetic
% problem 4 = conducting magnetic Remington shell casing
% problem 5 = Taufiq steel gun
% problem 6 = 72 anologue detenator
% problem 7 = belt buckle
% problem 8 = coin
% problem 9 = conducting magnetic torus
% problem 10 = Ansys Mesh (Sphere)

% Select which problem to run.
problem =5;
% Set polynomial degree of elements
pm =4;
% order of polynomial
maxorder = 10;
if pm > maxorder
    error('please run again with a lower order')
end


tic

% define all data for this problem
if problem==0
    % Simple plane wave propagation (for testing only)
    profun=@problem0;
    probdata = profun(pm);
elseif problem==1
    % unit conducing magnetic sphere, scaled by delta
    profun=@problem1;
    probdata = profun(pm);
elseif problem==2
    % unit spheroid, scaled by delta in a uniform H field
    profun=@problem1spheroid;
    probdata = profun(pm);
elseif problem==3
    % unit sphere in sphere sphere, scaled by delta
    profun=@problem1sins;
    probdata = profun(pm);
elseif problem==4
    % unit remington shell casing, scaled by delta
    profun=@problem1remington;
    probdata = profun(pm);
elseif problem==5
    % unit Taufiq gun, scaled by delta
    profun=@problem1taufiqgun;
    probdata = profun(pm);
elseif problem==6
    % unit landmine 72 analogue detanator, scaled by delta
    profun=@problem1landmine;
    probdata = profun(pm);
elseif problem==7
    % unit Belt Buckle, scaled by delta
    profun=@problem1beltbuckle;
    probdata = profun(pm);
elseif problem==8
    % unit 1 pound coin, scaled by delta
    profun=@problem1coin;
    probdata = profun(pm);
elseif problem==9
    % unit torus magnetic torus, scaled by delta
    profun=@problem1torus;
    probdata = profun(pm);
elseif problem==10
    % ANSYS Testing unit 1 sphere, scaled by delta
    profun=@problem1ansys;
    probdata = profun(pm);
end
% Job data
job = probdata.jb.job;
matflg = probdata.jb.matflg;
meshtype = probdata.jb.meshtype;
order = probdata.jb.order;
gorder = probdata.jb.gorder;
g1 =probdata.jb.g1;
shift = probdata.matr.shift;
delta=probdata.matr.delta;
nrhs=probdata.jb.nrhs;

% Compute maximum elemental matrix size
maxesize = 6*(maxorder+1) + 12*(maxorder-1) + 4*(maxorder-1)*(maxorder-2) +...
    2*(maxorder-1)*(maxorder-2) +(maxorder-1)*(maxorder-2)*(maxorder-3)/2;


if probdata.jb.meshtype==3
    meshdata = textread([job '.vol'], '%s','delimiter', '\n');
    disp('reading netgen mesh....')
    [nboun,bsido,nelem,intma,npoin,coord,cintma,mat,bcosurf,nsurf,quadlin]=meshinfo(meshdata,...
        meshtype,matflg,shift,delta);
    disp('...done')
else
      disp('reading ANSYS mesh....')
    [quadlin,intma,nelem,coord,mat,npoin,cintma]= ansysmeshinfo();
    disp('...done')
end


probdata.jb.quadlin=quadlin;
if quadlin==1
    % If linear geometry detected in mesh file fix linear interpolation of the
    % geometry
    gorder=0;
    probdata.jb.gorder=gorder;
end

% Flag those elements in the conductor
matc=zeros(nelem,1);
for i=1:nelem
    for j=probdata.matr.matcond
        if mat(i)==j
            matc(i)=1;
        end
    end
end

% assain the order to elements
if order > maxorder
    error(message('increae maxorder'));
else
    orderel = order*ones(nelem,1);
end

% size of local arrays exc interiors
esize = 6*(order+1);
if order >=2
    esize = esize +(12*(order-1))+(4*(order-1)*(order-2));
end
if esize > maxesize
    error(message('maxesize exceeded'));
end

% size of local arrays including interiors
esizet = esize;
if order>=3
    esizet=esizet+2*(order-1)*(order-2)+(order-2)*(order-1)*(order-3)/2 ;
end
if esizet > maxesize
    error(message('maxesize exceeded'));
end

% set all reference element types
[intma,cintma,eltype]=gettype(nelem,coord,intma,cintma,quadlin);
disp('Completed setting reference types')

% fix the centre of mesh to be in the centre of the object
[coord]=getcentreobj(coord,mat,matc,eltype,nelem,npoin,intma);
disp('Completed centering mesh')

% Automatic detection of boundaries
[bsido,bcosurf,nboun]=find_boundaries(npoin,nelem,intma,mat,matc);
disp('Completed automatic boundary detection')

% compute the numbering of the edges
[iside,nside,glob,help1,help2]=side(nelem,npoin,eltype,intma);
disp('Completed edge numbering')

% compute the numbering of the faces
[iface,help1,help2,help3,help4,globfa,cond,nface,bcedge]=faces_jin(npoin,...
    nelem,nboun,nside,help1,intma,bsido,glob,bcosurf);
disp('Completed face numbering')

% determine the unknown numbers
[unkz,unksid,unkfatp1,unkfatp2,unkfatp3,unkint,nef,npec,nunk,nunkt,...
    ordersid,orderfa,help1]=nounk(nelem,nside,nface,cond,glob,globfa,order,...
    orderel,bcedge,matc);

% zero edge coefficents
edgecof = zeros(nside,3*(gorder-1)+3);
% zero face coefficents
facecof = zeros(nface,(gorder*(gorder-1)/2-1)*3+3);

if gorder >=1
    % The currently available curved geometry types
    % 1 = not available
    % 2 = sphere exact geometry
    % 3 = not available
    % 4 = quadratic geometry file
    display(['The geometry type is ',num2str(g1)]);
    if g1 ==2
        [edgecof,facecof,pher]=getblendh1(nelem,nside,nface,intma,coord,...
            cond,glob,globfa,eltype,gorder,mat,probdata);
    elseif g1 ==4
        % This requires the mwah to be generated using quadractic option in
        % Netgen
        if quadlin==2
            [edgecof,facecof,pher]=getblendqh1(nelem,nside,nface,intma,coord,...
                cond,glob,globfa,eltype,gorder,mat,probdata,cintma,shift,delta,matc);
        else
            return
            error('using quadratic geometry')
        end
    end
    display('completed all bending');
end

display(['Program using p = ',num2str(order),'elements']);
display(['Size of elemental arrays = ',num2str(esizet)]);

% define nip = number of integration points in each direction
nipd = floor(((2*(order+1))+1+order)/2);
if nipd<=2
    nipd=2;
end

% Perform main matrix assembly
[rhs,stmtx,known,ephx1,ephy1,ephz1,ecphx1,ecphy1,ecphz1,ephx2,...
    ephy2,ephz2,ecphx2,ecphy2,ecphz2,intxi,inteta,intzeta,intw,nip,intfxi,...
    intfet,intfw,nipf,premtx,probdata,gphx1,gphy1,gphz1,gphx2,gphy2,gphz2,phh11,phh12]= eigenass(nelem,intma,coord,order,eltype,glob,unksid,...
    nside,cond,globfa,unkint,esizet,nipd,nunkt,nef,edgecof,facecof,...
    gorder,bcedge,npec,orderel,nrhs,unkfatp1,unkfatp2,unkfatp3,mat,unkz,probdata,matc);


% Solving Linear System
display('solving linear system')
if probdata.sol.regopt~=3
    [ Xzz,XEgEg,XFgFg, XFF, XIgIg, XII, nz, nheg, nhfg, nhf,nhig, nhi] = extract( unkz,unksid,unkfatp1,unkfatp2,unkfatp3,unkint,premtx,order,stmtx);
    tol=1e-7;
    maxit=1000;
    % Set up data for the preconditioner
    arg.Xzz=Xzz;
    arg.nz=nz;
    arg.XEgEg=XEgEg;
    arg.nheg=nheg;
    arg.XFgFg=XFgFg;
    arg.nhfg=nhfg;
    arg.XFF=XFF;
    arg.nhf=nhf;
    arg.XIgIg=XIgIg;
    arg.nhig=nhig;
    arg.XII=XII;
    arg.nhi=nhi;
end
if probdata.sol.regopt ==1
    % Solve linear system with iterative solve for application of the
    % action of the inverse of the gradient blocks
    % don't use eye as it generates a full matrix (uses lots of memory)
    I=[1:nunkt];
    J=[1:nunkt];
    X=ones(1,nunkt);
    preconlm=sparse(I,J,X,nunkt,nunkt);
    
    [L,U,P,Q] = lu(Xzz);
    disp('completed LU decomposition')
    arg.L=L;
    arg.U=U;
    arg.P=P;
    arg.Q=Q;
    if nhf > 0
        [Lf,Uf,Pf,Qf] = lu(XFF);
        disp('completed LU decomposition')
        arg.Lf=Lf;
        arg.Uf=Uf;
        arg.Pf=Pf;
        arg.Qf=Qf;
    end
    
    if nhi > 0
        [Li,Ui,Pi,Qi] = lu(XII);
        disp('completed LU decomposition')
        arg.Li=Li;
        arg.Ui=Ui;
        arg.Pi=Pi;
        arg.Qi=Qi;
    end
    sol=zeros(nunkt,nrhs);
    for i=1:nrhs
        [sol(:,i),flag,relres,iter,resvec] = gmres(stmtx,rhs(:,i),10,tol,maxit,preconlm,@(x)preconditioner(x,arg));
        if flag~=0
            disp('problem with gmres')
            flag
            semilogy(1:length(resvec),resvec/norm(rhs(:,i)),'-o');
            xlabel('Iteration number');
            ylabel('Relative residual');
        else
            disp('System solved succesfully');
            disp(['Number of iterations',num2str(length(resvec))]);
            semilogy(1:length(resvec),resvec/norm(rhs(:,i)),'-o');
            [[1:length(resvec)]',resvec/norm(rhs(:,i))]
            norm(rhs(:,i))
            xlabel('Iteration number');
            ylabel('Relative residual');
        end
    end
elseif probdata.sol.regopt==2
    % Solve linear system with direct LU solve for application of the
    % action of the inverse of the gradient blocks
    sol=zeros(nunkt,nrhs);
    I=[1:nunkt];
    J=[1:nunkt];
    X=ones(1,nunkt);
    preconlm=sparse(I,J,X,nunkt,nunkt);
    [L,U,P,Q] = lu(Xzz);
    disp('completed LU decomposition')
    arg.L=L;
    arg.U=U;
    arg.P=P;
    arg.Q=Q;
    if nheg > 0
        [Leg,Ueg,Peg,Qeg] = lu(XEgEg);
        disp('completed LU decomposition')
        arg.Leg=Leg;
        arg.Ueg=Ueg;
        arg.Peg=Peg;
        arg.Qeg=Qeg;
    end
    if nhfg > 0
        [Lfg,Ufg,Pfg,Qfg] = lu(XFgFg);
        disp('completed LU decomposition')
        arg.Lfg=Lfg;
        arg.Ufg=Ufg;
        arg.Pfg=Pfg;
        arg.Qfg=Qfg;
    end
    if nhf > 0
        [Lf,Uf,Pf,Qf] = lu(XFF);
        disp('completed LU decomposition')
        arg.Lf=Lf;
        arg.Uf=Uf;
        arg.Pf=Pf;
        arg.Qf=Qf;
    end
    if nhig>0
        [Lig,Uig,Pig,Qig] = lu(XIgIg);
        disp('completed LU decomposition')
        arg.Lig=Lig;
        arg.Uig=Uig;
        arg.Pig=Pig;
        arg.Qig=Qig;
    end
    
    if nhi > 0
        [Li,Ui,Pi,Qi] = lu(XII);
        disp('completed LU decomposition')
        arg.Li=Li;
        arg.Ui=Ui;
        arg.Pi=Pi;
        arg.Qi=Qi;
    end
    for i=1:nrhs
        [sol(:,i),flag,relres,iter,resvec] = gmres(stmtx,rhs(:,i),10,tol,maxit,preconlm,@(x)preconditioner_dir(x,arg));
        
        if flag~=0
            disp('problem with gmres')
            flag,resvec,relres
            semilogy(1:length(resvec),resvec/norm(rhs(:,i)),'-o');
            xlabel('Iteration number');
            ylabel('Relative residual');
        else
            disp('System solved succesfully');
            disp(['Number of iterations',num2str(length(resvec))]);
            semilogy(1:length(resvec),resvec/norm(rhs),'-o');
            xlabel('Iteration number');
            ylabel('Relative residual');
        end
    end
else
    % Direct Solve
    sol = stmtx\rhs;
end

% Finish up solve by adding Dirchlet DOFs to solution vector
[unkz,unksid,unkfatp1,unkfatp2,unkfatp3,unkint,sol]=streat(known,npec,sol,...
    nunkt,unksid,order,nside,nface,nelem,nrhs,unkfatp1,unkfatp2,unkfatp3,unkz,unkint);

if problem==0
    % this function is just for the purpose of checking the performance of the
    % hp-fem (simple wave propagation to confirm exponential convergence)
    hcurlerr =hcurl(nelem,gorder,eltype,edgecof,facecof,order,esizet,...
        intma,glob,globfa,unkz,unksid,unkfatp1,unkfatp2,unkfatp3,unkint,sol,coord,...
        nip,intxi,inteta,intzeta,intw,ephx1,ephy1,ephz1,ephx2,ephy2,ephz2,ecphx1,...
        ecphy1,ecphz1,ecphx2,ecphy2,ecphz2,probdata,gphx1,gphy1,gphz1,...
        gphx2,gphy2,gphz2,phh11,phh12)
else
    % Compute the magnetic rank 2 tensor
    ptensor =ptensorR2(nelem,gorder,eltype,edgecof,facecof,order,esizet,...
        intma,glob,globfa,sol,coord,unkz,unksid,unkfatp1,unkfatp2,unkfatp3,unkint,...
        nip,intxi,inteta,intzeta,intw,ephx1,ephy1,ephz1,ephx2,ephy2,ephz2,ecphx1,...
        ecphy1,ecphz1,ecphx2,ecphy2,ecphz2,probdata,matc,mat,gphx1,gphy1,gphz1,...
        gphx2,gphy2,gphz2,phh11,phh12);
    
    % In the case of a sphere compare with the analytical solution
    if probdata.mesh.erroroption==1
        m=exactspherev4(probdata);
        errorval=norm(ptensor-m*eye(3))/norm(m*eye(3));
        disp(['The relative error in the computed tensor=',num2str(errorval)])
    end
end

% Output the solution on a sub-divided mesh for purpose of visualisation
% using Paraview. Outputs .vtm and .vtu depending on job specified in
% problem file.
if probdata.mesh.vtkoption==1
    job=probdata.jb.job;
    pointvalue_multi(order,unkz,unksid,glob,globfa,eltype,...
        unkfatp1,unkfatp2,unkfatp3,unkint,esizet,sol,coord,intma,probdata,mat,gorder,edgecof,facecof,job);
end

toc

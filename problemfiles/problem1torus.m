function probdata=problem1(pm)

% Set the polynomial degree of the elements
% order = 0, 1, 2, 3
% nb order=0 refes to the standard linear/bilinear hat functions
% start adapting from this order

% include standard options (include BC, src definition etc)
probdata=[];
probdata=standoptions(probdata);

% Polynomial Degree Info Order
order=pm; 

% % plot the perturbed solution on a line
% % plotoption = 0 - do not plot
% % plotoption = 1 - plot
 plotoption = 0;

  
% Calculate the L2 norm of the error for the tensor
% % erroroption = 0 - do not calculate error
% % erroroption = 1 - calculate error
 erroroption = 0; % Exact E field not entered!
 
% % output the VTK file
% % vtkoption = 0 - do not output
% % vtkoption = 1 - output
 vtkoption =1;
 
 
probdata.mesh.plotoption=plotoption;
probdata.mesh.erroroption=erroroption;
probdata.mesh.vtkoption=vtkoption;

%job data------------------------------------------------------------------
job = 'torus'   % Job Name 
meshtype = 3;      % Mesh type 1= FLITE , 2=NG (old style), 3=NG (new style)
% % Blending Function Info
% gorder =0;        % order (if gorder > 0, quadlin =2 required)
% g1 = 4;            % type (g1=4 recomended)
% rin(1) = 1;        % raduis
% sufv(1)= 5;        % surfaceval
% gag = 2;           % Goagain  % 1- go again 2- once
% quadlin=0;         % quadratic geometry
% 
 probdata.jb.job = job;
 probdata.jb.meshtype = meshtype;
 probdata.jb.order=order;
% probdata.jb.gorder = gorder;
% probdata.jb.g1 = g1;
% probdata.jb.rin = rin(1);
% probdata.jb.sufv = sufv(1);
% probdata.jb.gag = gag(1);
% probdata.jb.nrhs=3;
% probdata.jb.quadlin = quadlin;

%material data-------------------------------------------------------------
nmat = 2;         % Number of materieals
muz = 1.256637061435917e-06; % Mu_z
epz = 0;                     % Ep_z
omega =133.5e2%              % Omega

% Mu_r, Ep_r, Sigma, J
% Mat 1
mu(1) = 1;
epl(1) = 0;
sigma(1) = 0.1;
jsrc(1,1:3) = [0  0  0];
% Mat 2
mu(2) = 1.5;
epl(2) = 0;
sigma(2) = 5.96e7;
jsrc(2,1:3) = [0  0  0];

% specify the material to be used a conductor conductors
% ie specify regions where gradients basis functions to be included
matcond=[2];

% In this case the mesh is for a unit sized object, it must be scaled (and
% repositioned)
delta = 0.01;  % Object size
shift=[0 0 0]; % Object shift

probdata.matr.nmat=nmat;
probdata.matr.muz=muz;
probdata.matr.epz=epz;
probdata.matr.omega=omega;
probdata.matr.mu=mu;
probdata.matr.epl=epl;
probdata.matr.sigma=sigma;
probdata.matr.jsrc=jsrc;
probdata.matr.delta=delta;
probdata.matr.shift=shift;
probdata.matr.matcond=matcond;
%-------------------------------------------------------------material data






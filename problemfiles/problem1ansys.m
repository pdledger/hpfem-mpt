function probdata=problem1ansys(pm)

% Set the polynomial degree of the elements
% order = 0, 1, 2, 3
% nb order=0 refes to edge elements with constant tangential coefficents

% include standard options (include BC, src definition etc)
probdata=[];
probdata=standoptions(probdata);

% Polynomial Degree Info Order
order=pm; 

% Output Options
%-------------------------------------------------------------------------
% % plot the perturbed solution on a line
% % plotoption = 0 - do not plot
% % plotoption = 1 - plot
 plotoption = 1;
 
% Calculate the L2 norm of the error for the tensor
% % erroroption = 0 - do not calculate error
% % erroroption = 1 - calculate error
 erroroption = 1; 
 
% % output the VTK file
% % vtkoption = 0 - do not output
% % vtkoption = 1 - output
 vtkoption =1;
 

probdata.mesh.plotoption=plotoption;
probdata.mesh.erroroption=erroroption;
probdata.mesh.vtkoption=vtkoption;

%-------------------------------------------------------------------------

%job data------------------------------------------------------------------
job = 'gsphere'   % Job Name 
meshtype = 4;      % Mesh type 1= FLITE , 2=NG (old style), 3=NG (new style)
                    %ansys

probdata.jb.job = job;
probdata.jb.meshtype = meshtype;
probdata.jb.order=order;
%------------------------------------------------------------------job data


probdata.sol.regopt =2 % 1- use regularisation with cg for gradent blocks
                       % 2- use regularisation with direct solve for gradient blocks
                       % 3- direct solve
%material data-------------------------------------------------------------
muz = 1.256637061435917e-06; % Mu_z
epz = 0;                     % Ep_z
omega =2*pi*133.5;           % Omega

% For each subdomain specify parameters
% Mu_r, Ep_r, Sigma, J
% Mat 1
mu(1) = 1;
epl(1) = 0;
sigma(1) = 0.1;
jsrc(1,1:3) = [0  0  0];

% Mat 2
mu(2) = 1.5
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

%-------------------------------------------------------------------------
% Blending Function Info (over write defaults)
% gorder =4;        % order (if gorder > 0, quadlin =2 required)
% g1 = 2;           % use exact geometry for a sphere
 rin(1) = 1;        % raduis
 sufv(1)= 5;        % surfaceval
 gag = 2;           % Goagain  % 1- go again 2- once
 svchk = 1;         % Check surface/volume as expecting a sphere
% 
% 
% probdata.jb.gorder = gorder;
% probdata.jb.g1 = g1;
 probdata.jb.rin = rin(1);
 probdata.jb.sufv = sufv(1);
 probdata.jb.gag = gag(1);
 probdata.mesh.svchk=svchk;

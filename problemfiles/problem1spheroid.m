function probdata=problem1spheroid(pm)

%-------------------------------------------------------------------------
probdata=[];
% include standard options (include BC, src definition etc)
probdata=standoptions(probdata);

% Polynomial Degree Info Order
order=pm; 

% % plot the perturbed solution on a line
% % plotoption = 0 - do not plot
% % plotoption = 1 - plot
plotoption = 1;
 

% Calculate the L2 norm of the error for the tensor
% % erroroption = 0 - do not calculate error
% % erroroption = 1 - calculate error
erroroption = 0; 
 
% % output the VTK file
% % vtkoption = 0 - do not output
% % vtkoption = 1 - output
vtkoption =1;
 
probdata.mesh.plotoption=plotoption;
probdata.mesh.erroroption=erroroption;
probdata.mesh.vtkoption=vtkoption;

%job data------------------------------------------------------------------
job = 'ellipsemod'   % Job Name 
meshtype = 3;        % Mesh type 1= FLITE , 2=NG (old style), 3=NG (new style)

probdata.jb.job = job;
probdata.jb.meshtype = meshtype;
probdata.jb.order=order;
% %------------------------------------------------------------------job data


%material data-------------------------------------------------------------
nmat = 2;         % Number of materieals
muz = 1.256637061435917e-06; % Mu_z
epz = 0;                     % Ep_z
omega =133.5;              % Omega

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

%-------------------------------------------------------------------------
% Blending Function Info (overwrite defaults)
gorder =0;        % order (if gorder > 0, quadlin =2 required)
probdata.jb.gorder = gorder;

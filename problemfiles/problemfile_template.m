function probdata=problemfile_template(pm)

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
 erroroption = 0; 
 
% % output the VTK file
% % vtkoption = 0 - do not output
% % vtkoption = 1 - output
 vtkoption =1;
 
probdata.mesh.plotoption=plotoption;
probdata.mesh.erroroption=erroroption;
probdata.mesh.vtkoption=vtkoption;

%-------------------------------------------------------------------------

%job data------------------------------------------------------------------
% Input needed here
job = '%insert job name here%'   % Job Name
% eg  job='sphere11' if the mesh file is sphere11.vol
meshtype = 3;                    % Mesh type 1= FLITE , 2=NG (old style), 3=NG (new style)
% no change to mesh type required if using NETGEN.


probdata.jb.job = job;
probdata.jb.meshtype = meshtype;
probdata.jb.order=order;
%------------------------------------------------------------------job data

%material data-------------------------------------------------------------
muz = 1.256637061435917e-06; % Mu_z
epz = 0;                     % Ep_z

%Input needed here,
omega = %insert frequency in rad/s;           
% eg if omega is 133.5 rad/s enter 
% omega =133.5


% Input needed here
% For each subdomain specify parameters
% Mu_r, Ep_r, Sigma, J
% Mat 1

%eg if the first subdomain is non-conducting enter
%mu(1) = 1;
%epl(1) = 0;
%sigma(1) = 0.1;
%jsrc(1,1:3) = [0  0  0];

% and the second subdomain is the conducting object with mu_r =1.5 and
% sigma =5.96e7 S/m enter
%mu(2) = 1.5
%epl(2) = 0;
%sigma(2) = 5.96e7;
%jsrc(2,1:3) = [0  0  0];


% Input needed here
% specify the material to be used a conductor conductors
% ie specify regions where gradients basis functions to be included
matcond= % insert an array contain the list of conductings subdomains.
% eg if the only conducting subdomain is the second one enter
% matcond=[2];

% Input needed here
% In this case the mesh is for a unit sized object, it must be scaled (and
% repositioned)
delta = % enter the object scaling
% specify what scaling to apply to go from the O(1) object specified in the
% mesh to the actual object. Eg if a conducting object of O(1) is specified
% in the mesh but the true object has dimension of O(0.01) enter
%delta = 0.01;  % Object size
shift=[0 0 0]; % Object shift

%probdata.matr.nmat=nmat;
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


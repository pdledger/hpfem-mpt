% Function specifying the problem details for a Remington 222 shell casing

function probdata=problem1remington(pm)

% include standard options (include BC, src definition etc)
probdata=[];
probdata=standoptions(probdata);

% Polynomial Degree Info Order
order=pm; 

% % plot the perturbed solution on a line
% % plotoption = 0 - do not plot
% % plotoption = 1 - plot
 plotoption = 1;
 
 
% l2 norm norm of the tensor
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
job = 'remington222zpaper'   % Job Name 
meshtype = 3;      % Mesh type 1= FLITE , 2=NG (old style), 3=NG (new style)

probdata.jb.job = job;
probdata.jb.meshtype = meshtype;
probdata.jb.order=order;

%material data-------------------------------------------------------------
muz = 1.256637061435917e-06; % Mu_z
epz = 0;                     % Ep_z
omega =133.5e2;              % Omega

% Mu_r, Ep_r, Sigma, J
% Free space
for i=1:4
mu(i) = 1;
epl(i) = 0;
sigma(i) = 10;%0.1; % thin conductor requires increased regularisation paramter
jsrc(i,1:3) = [0  0  0];
end
% Object
for i=5:8
mu(i) = 1.0;
epl(i) = 0;
sigma(i) = 1.5e7;
jsrc(i,1:3) = [0  0  0];
end
% specify the material to be used a conductor conductors
% ie specify regions where gradients basis functions to be included
matcond=[5:8];

% In this case the mesh is for a unit sized object, it must be scaled (and
% repositioned)
delta = 0.001;  % Object size
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
%-------------------------------------------------------------material data
% Blending Function Info (over write defualts)
gorder =0;        % order (if gorder > 0, quadlin =2 required)
probdata.jb.gorder = gorder;


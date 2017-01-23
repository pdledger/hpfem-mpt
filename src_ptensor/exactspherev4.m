% Function to compute the exact polarizability tensor for a conducting sphere

function m=exactspherev4(probdata)

matcond=probdata.matr.matcond;


% free space permeability (in H/m)
mu0=probdata.matr.muz;

% Angular Frequency (in rad/s)
omega=probdata.matr.omega;

% object relative permeability
muv=probdata.matr.mu;
mu=muv(matcond);

% object conductivity
sigmav=probdata.matr.sigma;
sigma=sigmav(matcond);

alpha=probdata.matr.delta;

epsilon=0;

%find the polarisation tensor
m=polarisation(omega,epsilon,sigma,mu*mu0,mu0,alpha);
disp(['The exact polarisation tensor is diagonal with coefficent',num2str(m)])


function m=polarisation(omega,epsilon,sigma,mu,mu0,alpha)

r0=alpha;

% compute terms
k=sqrt(omega^2*epsilon*sigma  +i*mu*sigma*omega);

js0kr=sqrt(pi/(2*k*r0))*besselj(1/2,k*r0);
js1kr=sqrt(pi/(2*k*r0))*besselj(1+1/2,k*r0);
js2kr=sqrt(pi/(2*k*r0))*besselj(2+1/2,k*r0);

% object's polarisation tensor
m=2*pi*alpha^3* (2*(mu-mu0)*js0kr+(2*mu+mu0)*js2kr)/( (mu+2*mu0)*js0kr + (mu-mu0)*js2kr);

function ptf=perturbfield(x,z,m,H0);

r=x-z;
rad=norm(r);
rhat=r/rad;

D2G=1/4/pi/rad^3*(3*rhat*rhat'-eye(3));

ptf=m*D2G*H0;


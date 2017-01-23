% Set up Standard Options for BC's, src data etc for an MPT calculation
function probdata=standoptions(probdata)


%-------------------------------------------------------------------------
% Not used
% % plot the  voltage  on a line
% % voltoption = 0 - do not plot
% % voltoption = 1 - plot
 voltoption = 0; 

% % plot the  H_alpha-H_0 - D2G PT H_0  on a line
% % darkooption = 0 - do not plot
% % darkooption = 1 - plot
 darkooption = 0; 

% % plot the  voltage  on a line
% % voltoption = 0 - do not plot
% % voltoption = 1 - plot
 voltoption = 0; 
 
probdata.mesh.voltoption=voltoption;
probdata.mesh.darkooption=darkooption;
 
% Not used.
% %-----------------------------------------------------------------------
%  %interpolation points for perturbed H
%  
% % %interpolation number
%  N = 50;
%  probdata.mesh.N = N;
%  
% 
% starp = [0.0,0,0.01];   % startpoin of the line
%  overp =[0.0,0,0.08];   % endpoint of the line
%       point = starp;   
%      deltax = (overp(1)-starp(1))/(N-1);
%      deltay = (overp(2)-starp(2))/(N-1);
%      deltaz = (overp(3)-starp(3))/(N-1);
% 
%      for i = 2:N
%          point = [point;starp(1)+deltax*(i-1) starp(2)+deltay*(i-1) starp(3)+deltaz*(i-1)];    %coordinate of interpolation points
%      end 
%      
%  probdata.mesh.point = point;
% %-----------------------------------------------------
% % define coil positions (centres) for measuring the voltage
% 
% Npts=10;
% dphi=2*pi/Npts;
% coilpos=[];
% % coil aligned with z axis
% coilrad=0.01;
% coilcenz=[0.01 0.02 0.03 0.04 0.05];
% Ncoils=length(coilcenz);
% for r=coilcenz
%     for i = 1:Npts
%          coilpos = [coilpos; coilrad*cos((i-1)*dphi) coilrad*sin((i-1)*dphi) r];    %coordinate of coils
%     end 
% end
% 
% probdata.coils.Npts=Npts;       % number of points defining the coil
% probdata.coils.Ncoils=Ncoils;     % number of coils
% probdata.coils.coilpos=coilpos; % locations where field is to be computed
% probdata.coils.coilrad=coilrad; % coil radii
% probdata.coils.coilcenz=coilcenz; % coil centre



% solver options ---------------------------------------------------------

probdata.sol.regopt =2 % 1- use regularisation with cg for gradent blocks
                       % 2- use regularisation with direct solve for gradient blocks
                       % 3- direct solve

                       
probdata.jb.nrhs=3;  % 3 right hand sides
%--------------------------------------------------------------------------
% Standard mesh options
matflg = 1;        % Mat Flag from file
ppus   = 1 ;       % Parrellel Processors (not used)
% Polynomial Degree Info
ptype = 2;         % (Type)

probdata.jb.matflg=matflg;
probdata.jb.ppus = ppus;
probdata.jb.ptype = ptype; 

%-------------------------------------------------------------------------
% Standard Geometry Options

% Blending Function Info
gorder =2;        % order (if gorder > 0, quadlin =2 required)
g1 = 4;            % type (g1=4 recomended)
sufv(1)= 5;        % surfaceval
gag = 2;           % Goagain  % 1- go again 2- once
svchk = 0;         % Do check surface area /volume as not expecting a sphere

probdata.jb.gorder = gorder;
probdata.jb.g1 = g1;
probdata.jb.sufv = sufv(1);
probdata.jb.gag = gag(1);
probdata.mesh.svchk=svchk;
%-------------------------------------------------------------------------

% define boundary data

% current bctype
% bctype 5 Object (interface between object and free space)
% bctype 3 Neumann type (expected for this problem)
% bctype 2 Dirichlet (none here)

% Note that boundaries will be automatically detected so don't need to be
% specified in the mesh file.

% define the function handle for Dirichlet BC's
probdata.es.dirfunarg=[];
probdata.es.dirfun=@esproblemdir;

% define the source terms
probdata.es.srcfunarg=[];
probdata.es.srcfun=@esproblemsrc;

% define the function handle for Neumann BC's
H0=[0;0;1];
arg.H0=H0;
muz = 1.256637061435917e-06; % Mu_z
arg.muz=muz;
probdata.es.neufun=@esproblemneu;
probdata.es.neufunarg=arg;

% exact solution
%arg.mtensor=mtensor
% probdata.es.exactfun=@esproblemexact; % not entered
% probdata.es.exactfunarg=arg;
% 
% probdata.es.exactcurlfun=@esproblemexactcurl;
% probdata.es.exactcurlfunarg=arg;

%--------------------------------------------------------------------------

% Dirichlet Boundary Conditio----------------------------------------------
function e=esproblemdir(x,y,z,index,arg)       % DBC for es (4*pi*10e-7)^(0.5)*

domega = arg;
if index==2
    % this would be non-zero for non zero type BC's
    e=zeros(3,3);
else
    e=zeros(3,3);
end

% Neumann Boundary Conditio------------------------------------------------

function curle=esproblemneu(x,y,z,index,arg)

if index==5
    curle=-2*eye(3);
else
    curle=zeros(3,3);
end

% Src term-----------------------------------------------------------------

function src=esproblemsrc(x,y,z,dum,arg)               % Source term for es
kappa=arg.kappa;
xyz=[x;y;z];
IT=eye(3);
for i=1:3
    uu=IT(:,i);
    out=-kappa*cross(uu,xyz);
    src(:,i)=out;
end

% % Analytical Solution Field------------------------------------------------
% % not given for this problem
% function [out]=esproblemexact(x,y,z,arg)           
% out=zeros(3,1);
% 
% 
% % compute perturbed field
% function ptf=perturbfield(xv,zv,mtensor,H0);
% r=xv-zv;
% rad=norm(r);
% rhat=r/rad;
% 
% D2G=1/4/pi/rad^3*(3*rhat*rhat'-eye(3));
% 
% ptf=D2G*(mtensor*H0);
%  
% % Analytical Solution for B=curl A =mu H-----------------------------------
% % Gives total field as output
% function [curle,H0]=esproblemexactcurl(x,y,z,arg)         
% 
% muz = arg.muz;
% H0=arg.H0;
% mtensor=arg.mtensor;
% rpos=[x,y,z]';
% zpos=[0,0,0]';
% % perturbed magnetic field
% ptf=perturbfield(rpos,zpos,mtensor,H0);
% % back ground field constant so not recomputed
% curle=(H0+ptf)*muz;

%--------------------------------------------------------------------------
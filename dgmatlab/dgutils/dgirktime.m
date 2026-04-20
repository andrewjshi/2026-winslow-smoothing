function u=dgirktime(u,msh,data,phys,jacname,tpars,spars,solver)
%DGIRKTIME  DG Implicit Runge-Kutta Time Stepping
%
%    Syntax: u=dgirktime(u,msh,data,phys,jacname,tpars,spars,solver)
%
%    u:                Initial solution (input)
%    u:                Final solution (output)
%    msh, data, phys:  DG problem specification
%    jacname:          DG Residual/Jacobian function
%    tpars:            = [dt, nsteps, nstages]
%                dt:     Timestep
%            nsteps:     Number of timesteps
%           nstages:     Number of implicit Runge-Kutta stages
%                        Order of accuracy = nstages + 1
%    spars:            = [newtontol, newtonmaxiter, linsystol,
%                         linsysmaxiter, reportlevel]
%         newtontol:     Residual tolerance in Newton iterations (1e-8)
%     newtonmaxiter:     Maximum number of Newton steps           (20)
%         linsystol:     Residual tolerance in iterative solver  (1e-5)
%     linsysmaxiter:     Maximum number of Krylov iterations     (1000)
%       reportlevel:     Level of diagnostic output in linear solver (0)
%    solver:           Solver command
%                      = 'direct'     : Direct solver (default)
%                        'gmres j'    : GMRES, block Jacobi precond.
%                        'gmres i'    : GMRES, block ILU(0) precond.
%                        'gmres 0dpi' : GMRES, ILU(0)-p0 precond.
%                        'gmres'      :       -    "    -

if nargin<7 | isempty(spars), spars=[]; end
if nargin<8 | isempty(solver), solver='direct'; end
[spars,solvecommand,precond]=parsesolver(spars,solver);

dt=tpars(1); nsteps=tpars(2); nstages=tpars(3);
newtontol=spars(1); maxiter=spars(2);

switch nstages
 case 1
%  A=1/2; b=1;
   A=1; b=1;
 case 2
%  A=[(1+1/sqrt(3))/2, 0; -1/sqrt(3), (1+1/sqrt(3))/2]; b=[1;1]/2;
   al=1-1/sqrt(2);
   A=[al,0; 1-al,al];
   b=[1-al,al];
 case 3
%  al=2*cos(pi/18)/sqrt(3);
%  A=(1+al)/2*eye(3)+[0,0,0;-al/2,0,0;1+al,-(1+2*al),0];
%  b=[1/6/al^2,1-1/3/al^2,1/6/al^2];
   al=0.435866521508459;
   t2=(1+al)/2;
   b1=-(6*al^2-16*al+1)/4;
   b2=(6*al^2-20*al+5)/4;
   A=[al,0,0; t2-al,al,0; b1,b2,al];
   b=[b1,b2,al];
end
c=sum(A,2);

r=feval(jacname,u,msh,data,phys);
N=size(r,2);
if ~isfield(phys,'time'), phys.time=0.0; end
M=mkmassmatrix(msh,data,N,true);
for ii=1:nsteps
  time0=phys.time;
  k=zeros([size(u,1),N,size(u,3),nstages]);
  fprintf('Timestep %4d\n',ii);
  for is=1:nstages
    phys.time=time0+dt*c(is);
    fprintf('  Stage %d\n',is);
    if A(is,is)==0
      cu=u;
      cu(:,1:N,:)=cu(:,1:N,:)+dt*kdot(k,A(is,1:is-1));
      k(:,:,:,is)=dgmassinv(feval(jacname,cu,msh,data,phys),msh,data,0);
    else
      cu=u;
      cu(:,1:N,:)=cu(:,1:N,:)+dt*kdot(k,A(is,1:is));
      R=feval(jacname,cu,msh,data,phys);
      fprintf('  Newton iter %2d, Res = %18.14f',0,max(abs(R(:))));
      for jj=1:maxiter
        clear DA OA
        cu=u;
        cu(:,1:N,:)=cu(:,1:N,:)+dt*kdot(k,A(is,1:is));
        [R,DA,OA]=feval(jacname,cu,msh,data,phys);
        R=dgmassinv(k(:,:,:,is),msh,data,1)-R;
        DA=M-dt*A(is,is)*DA;
        OA=-dt*A(is,is)*OA;
        
        [dk,outpars]=dglinsys(msh,data,DA,OA,R, ...
                              precond,spars(3:5),solvecommand);
        if strcmp(solver,'direct'), fprintf('\n'); else
          fprintf(', Krylov iters = %4d, res = %5.2g\n',outpars(2),outpars(1));
        end
        if any(isnan(dk(:))), error('NaN in solution.'); end
        k(:,:,:,is)=k(:,:,:,is)-dk;

        clear DA OA
        cu=u;
        cu(:,1:N,:)=cu(:,1:N,:)+dt*kdot(k,A(is,1:is));
        R=feval(jacname,cu,msh,data,phys);
        R=dgmassinv(k(:,:,:,is),msh,data,1)-R;
        normR=max(abs(R(:)));
        fprintf('  Newton iter %2d, Res = %18.14f',jj,normR);
        if normR<newtontol, break; end
      end
      fprintf('\n');
      if normR>=newtontol
        warning('No convergence.');
      end
    end
  end
  
  u(:,1:N,:)=u(:,1:N,:)+dt*kdot(k,b);
  phys.time=time0+ii*dt;
end

function v=kdot(k,a)

v=0*k(:,:,:,1);
for i=1:numel(a)
  v=v+a(i)*k(:,:,:,i);
end

function [spars,solvecommand,precond]=parsesolver(spars,solver)

if length(spars)<1, spars=[spars,1e-8]; end
if length(spars)<2, spars=[spars,20]; end
if length(spars)<3, spars=[spars,1e-5]; end
if length(spars)<4, spars=[spars,1000]; end
if length(spars)<5, spars=[spars,0]; end

ix=find(solver==' ');
if isempty(ix)
  solvecommand=solver;
  precond='';
else
  solvecommand=solver(1:ix-1);
  precond=solver(ix+1:end);
end

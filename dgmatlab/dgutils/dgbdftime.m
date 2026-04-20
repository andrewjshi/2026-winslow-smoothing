function u=dgbdftime(u,msh,data,phys,jacname,tpars,spars,solver,pltpars,filename)
%DGBDFTIME  BDF Time Stepping
%
%    Solve nonlinear problem with damped Newton iterations.
%
%    Syntax: u=dgbdftime(u,msh,data,phys,jacname,tpars,spars,solver,pltpars,filename)
%
%    u:                Initial guess (input)
%    u:                Solution (output)
%    msh, data, phys:  DG problem
%    jacname:          Jacobian/Residual function
%    tpars:            = [dt, nsteps, pltfreq]
%                dt:     Timestep
%            nsteps:     Number of timesteps
%           pltfreq:     Plotting frequency (inf for no plotting)
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
%
%    Optional plotting parameters:
%
%    pltpars:          Parameters to dgplot (cell array)
%                      Empty for no plotting
%    filename:         Filename to save solution (for recovery)

if nargin<7 | isempty(spars), spars=[]; end
if nargin<8 | isempty(solver), solver='direct'; end
if nargin<9, pltpars={}; end
if nargin<10, filename=''; end

if length(spars)<1, spars=[spars,1e-8]; end
if length(spars)<2, spars=[spars,20]; end
if length(spars)<3, spars=[spars,1e-5]; end
if length(spars)<4, spars=[spars,1000]; end
if length(spars)<5, spars=[spars,0]; end

if length(tpars)<3, tpars(3)=inf; end

ix=find(solver==' ');
if isempty(ix)
  solvecommand=solver;
  precond='';
else
  solvecommand=solver(1:ix-1);
  precond=solver(ix+1:end);
end

dt=tpars(1);
nsteps=tpars(2);
pltfreq=tpars(3);

newtontol=spars(1);
maxiter=spars(2);

dim=size(msh.s,2)-1;

if ~isfield(phys,'time'), phys.time=0.0; end

for ii=1:nsteps
  fprintf('Timestep %4d\n',ii);
  phys.time=phys.time+dt;

  clear DA OA
  R=feval(jacname,u,msh,data,phys);
  normR=dt*max(abs(R(:)));
  if ii==1, N=size(R,2); M=mkmassmatrix(msh,data,N,true); end

  u0=u;
  fprintf('  Newton iter %2d, Res = %18.14f',0,normR);
  for jj=1:maxiter
    clear DA OA
    [R,DA,OA]=feval(jacname,u,msh,data,phys);
    R=dgmassinv(u(:,1:N,:)-u0(:,1:N,:),msh,data,1)-dt*R;
    DA=M-dt*DA;
    OA=-dt*OA;
   
    [du,outpars]=dglinsys(msh,data,DA,OA,R,precond,spars(3:5),solvecommand);
    if strcmp(solver,'direct')
      fprintf('\n');
    else
      fprintf(', Krylov iters = %4d, res = %5.2g\n',outpars(2),outpars(1));
    end
    if any(isnan(du(:))), error('NaN in solution.'); end
    u(:,1:N,:)=u(:,1:N,:)-du;
    
    clear DA OA
    R=feval(jacname,u,msh,data,phys);
    R=dgmassinv(u(:,1:N,:)-u0(:,1:N,:),msh,data,1)-dt*R;
    normR=max(abs(R(:)));
    fprintf('  Newton iter %2d, Res = %18.14f',jj,normR);
    if normR<newtontol, break; end
  end
  fprintf('\n');
  if normR>=newtontol
    warning('No convergence.');
  end
  
  if ~isempty(pltpars) & mod(ii,pltfreq)==0
    dgplot(msh,u,pltpars{:});
  end
  if ~isempty(filename)
    save(filename,'u');
  end
end

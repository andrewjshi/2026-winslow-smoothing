function u=dgsolve(u,msh,data,phys,jacname,spars,solver,pltpars)
%DGSOLVE  Damped Newton solver
%
%    Solve nonlinear problem with damped Newton iterations.
%
%    Syntax: u=dgsolve(u,msh,data,phys,jacname,spars,solver,pltpars)
%
%    u:                Initial guess (input)
%    u:                Solution (output)
%    msh, data, phys:  DG problem
%    jacname:          Jacobian/Residual function
%    spars:            = [newtontol, newtonmaxiter, linsystol,
%                         reportlevel, damped]
%         newtontol:     Residual tolerance in Newton iterations (1e-8)
%     newtonmaxiter:     Maximum number of Newton steps           (20)
%         linsystol:     Residual tolerance in iterative solver  (1e-5)
%     linsysmaxiter:     Maximum number of Krylov iterations     (1000)
%    linreportlevel:     Level of diagnostic output in linear solver (0)
% newtonreportlevel:     Level of diagnostic output in Newton solver (1)
%            damped:     Use damped Newton solver (0)
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

if nargin<6 | isempty(spars), spars=[]; end
if nargin<7 | isempty(solver), solver='direct'; end
if nargin<8, pltpars={}; end

if length(spars)<1, spars=[spars,1e-8]; end
if length(spars)<2, spars=[spars,20]; end
if length(spars)<3, spars=[spars,1e-5]; end
if length(spars)<4, spars=[spars,1000]; end
if length(spars)<5, spars=[spars,0]; end
if length(spars)<6, spars=[spars,1]; end
if length(spars)<7, spars=[spars,0]; end

ix=find(solver==' ');
if isempty(ix)
  solvecommand=solver;
  precond='';
else
  solvecommand=solver(1:ix-1);
  precond=solver(ix+1:end);
end

newtontol=spars(1);
maxiter=spars(2);
damped=spars(7);
rep=spars(6);

dim=size(msh.s,2)-1;

R=feval(jacname,u,msh,data,phys);
N=size(R,2);

normR=max(abs(R(:)));
if rep>=1, fprintf('Iter %2d, Res = %18.14f',0,normR); end
if normR<newtontol
  if rep>=1, fprintf('\n'); end
  return;
end

for ii=1:maxiter
  % Solve
  clear DA OA
  [R,DA,OA]=feval(jacname,u,msh,data,phys);
  [du,outpars]=dglinsys(msh,data,DA,OA,R,precond,spars(3:5),solvecommand);
  if rep>=1
    if strcmp(solver,'direct')
      fprintf('\n');
    else
      fprintf(', Krylov iters = %4d, res = %5.2g\n',outpars(2),outpars(1));
    end
  end
  if any(isnan(du(:))), error('NaN in solution.'); end

  if damped
    % Find stepsize for damped Newton
    alpha=1;
    u0=u;
    while 1
      u(:,1:N,:)=u0(:,1:N,:)-alpha*reshape(du,size(R));
      R=feval(jacname,u,msh,data,phys);
      if max(abs(R(:)))>normR
        alpha=alpha/2;
        if rep>=1
          fprintf('  alpha = %f\n',alpha);
        end
        if alpha<1e-3, error('Too small alpha.'); end
      else
        break;
      end
    end
  else
    u(:,1:N,:)=u(:,1:N,:)-reshape(du,size(R));
  end
  
  if ~isempty(pltpars)
    dgplot(msh,u,pltpars{:});
  end
  
  % Compute residual
  R=feval(jacname,u,msh,data,phys);
  normR=max(abs(R(:)));
  if rep>=1
    fprintf('Iter %2d, Res = %18.14f',ii,normR);
  end
  if normR<newtontol, break; end
end
if rep>=1
  fprintf('\n');
end

if normR>=newtontol
  warning('No convergence.');
end

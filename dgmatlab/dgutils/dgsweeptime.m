function [u,stat]=dgsweeptime(u,msh,data,phys,jacname,tpars,spars,solver,pltpars,filename)
%DGSWEEPTIME  BDF Time Stepping, sweeping timestep
%
%    Implicit timestepping using the BDF method (first order) with
%    timesteps varying from dt1 to dt2. If convergence is "fast" (<=it1
%    iterations), the timestep is multiplied by fac. If convergence is
%    "normal" (<=it2 iterations), the timestep is unchanged. If no
%    convergence in it2 iterations, solution is discarded and timestep is
%    divided by fac.
%
%    Syntax: u=dgsweeptime(u,msh,data,phys,jacname,tpars,spars,solver,pltpars)
%
%    tpars:            = [dt1, dt2, fac, it1, it2, mindt]
%    spars:            = [newtontol, newtonmaxiter, linsystol,
%                         linsysmaxiter, reportlevel, maxsteps]
%      Like DGBDFTIME, but actual Newton tolerance is newtontol*dt, default
%      newtontol=1e-5, and newtonmaxiter is ignored (controlled by it1,it2).
%
%    See DGBDFTIME for other parameters.

if nargin<7 | isempty(spars), spars=[]; end
if nargin<8 | isempty(solver), solver='direct'; end
if nargin<9, pltpars={}; end
if nargin<10, filename=''; end

if length(tpars)<6, tpars(6)=tpars(1)*1e-1; end

if length(spars)<1, spars=[spars,1e-5]; end
if length(spars)<2, spars=[spars,20]; end
if length(spars)<3, spars=[spars,1e-5]; end
if length(spars)<4, spars=[spars,1000]; end
if length(spars)<5, spars=[spars,0]; end
if length(spars)<6, spars=[spars, 1000]; end

ix=find(solver==' ');
if isempty(ix)
  solvecommand=solver;
  precond='';
else
  solvecommand=solver(1:ix-1);
  precond=solver(ix+1:end);
end

dt1=tpars(1);
dt2=tpars(2);
fac=tpars(3);
it1=tpars(4);
it2=tpars(5);
mindt=tpars(6);

newtontol=spars(1);
maxsteps=spars(6);

dim=size(msh.s,2)-1;

fprintf('\n--- DGSWEEPTIME ---\n');
fprintf('dt1=%g, dt2=%g, fac=%g, it1=%d, it2=%d\n', ...
        dt1,dt2,fac,it1,it2);

stat.newtoniters=[];
stat.dts=[];
stat.acceptedsteps=[];
stat.T=0;

ii=0;
dt=dt1;
tic;
while 1
  ii=ii+1;
  fprintf('\n - Step %d, dt=%g\n',ii,dt);

  clear DA OA
  R=feval(jacname,u,msh,data,phys);
  normR=dt*max(abs(R(:)));
  if ii==1, N=size(R,2); M=mkmassmatrix(msh,data,N,true); end

  stat.dts(ii)=dt;
  u0=u;
  fprintf('  Newton iter %2d, Res = %18.14f',0,normR);
  normR0=normR;
  for jj=1:it2
    clear DA OA
    [R,DA,OA]=feval(jacname,u,msh,data,phys);
    R=dgmassinv(u(:,1:N,:)-u0(:,1:N,:),msh,data,1)-dt*R;
    DA=M-dt*DA;
    OA=-dt*OA;

    try
      [du,outpars]=dglinsys(msh,data,DA,OA,R,precond,spars(3:5),solvecommand);
    catch
      du=nan;
    end
    if strcmp(solver,'direct')
      fprintf('\n');
    else
      fprintf(', Krylov iters = %4d, res = %5.2g\n',outpars(2),outpars(1));
    end
    if any(isnan(du(:))), normR=nan; break; end
    u(:,1:N,:)=u(:,1:N,:)-du;
    
    clear DA OA
    R=feval(jacname,u,msh,data,phys);
    R=dgmassinv(u(:,1:N,:)-u0(:,1:N,:),msh,data,1)-dt*R;
    normR=max(abs(R(:)));
    fprintf('  Newton iter %2d, Res = %18.14f',jj,normR);
    if normR<newtontol*dt, break; end
    if normR>100*normR0, break; end
  end
  fprintf('\n');
  if isnan(normR) || normR>100*normR0 | normR>=newtontol*dt
    fprintf('No convergence: Discarding u, decreasing dt\n');
    stat.acceptedsteps(ii)=false;
    u=u0;
    dt=dt/fac;
  else
    stat.T=stat.T+dt;
    stat.acceptedsteps(ii)=true;
    if jj<=it1
      fprintf('Fast convergence: Keeping u, increasing dt\n');
      dt=dt*fac;
    else
      fprintf('Normal convergence: Keeping u, keeping dt\n');
    end
  end
  stat.newtoniters(ii)=jj;

  if ~isempty(pltpars)
    dgplot(msh,u,pltpars{:});
  end
  if ~isempty(filename)
    save(filename,'u');
  end

  if dt>dt2
    fprintf('\n--- DGSWEEPTIME FINISHED SUCCESSFULLY ---\n');
    break;
  end
  if ii>=maxsteps
    fprintf('\n--- DGSWEEPTIME FINISHED WITHOUT REACHING FINAL DT ---\n');
    break;
  end
  if dt<mindt
    fprintf('\n--- DGSWEEPTIME TERMINATED, TOO SMALL DT ---\n');
    break;
  end
end
stat.executiontime=toc;

fprintf('Execution time: %gs\n',stat.executiontime);
fprintf('Number of attempted steps: %d\n',length(stat.dts));
fprintf('Number of accepted steps: %d\n',sum(stat.acceptedsteps));
fprintf('Number of Newton iterations: %d\n',sum(stat.newtoniters));
fprintf('Total integrated time: %g\n',stat.T);
fprintf('\n');

function u=dgrktime(u,msh,data,phys,resname,dt,ipars,pltpars,moviefn)
%DGRKTIME  DG Runge-Kutta Timestepping
%
%    Take nsteps with RK4 explicit timestepping scheme.
%
%    Syntax: u=dgrktime(u,msh,data,phys,resname,dt,ipars,pltpars,moviefn)
%
%    u:                Initial solution (input)
%    u:                Final solution (output)
%    msh, data, phys:  DG problem
%    resname:          Residual function
%    dt:               Stepsize
%    ipars:            = [nsteps,pltfreq]
%            nsteps:     Total number of timesteps
%           pltfreq:     Plotting frequency (inf for no plotting)
%
%    Optional plotting parameters:
%
%    pltpars:          Parameters to dgplot (cell array)
%    moviefn:          Directory and partial filename for images

nsteps=ipars(1);
pltfreq=ipars(2);
dim=size(msh.s,2)-1;
if nargin<8, pltpars={}; end
if nargin<9, moviefn=''; end

if ~isfield(phys,'time')
  phys.time=0.0;
end
for ii=1:nsteps
  fprintf('Timestep %4d\n',ii);

  k1=dt*dgmassinv(feval(resname,u      ,msh,data,phys),msh,data);
  phys.time=phys.time+dt/2;
  k2=dt*dgmassinv(feval(resname,u+k1/2 ,msh,data,phys),msh,data);
  k3=dt*dgmassinv(feval(resname,u+k2/2 ,msh,data,phys),msh,data);
  phys.time=phys.time+dt/2;
  k4=dt*dgmassinv(feval(resname,u+k3   ,msh,data,phys),msh,data);
  u=u+(k1+2*k2+2*k3+k4)/6;
  if any(isnan(u(:))), error('NaN in solution.'); end

  if pltfreq & mod(ii,pltfreq)==0
    dgplot(msh,u,pltpars{:});
    if ~isempty(moviefn)
      dgimgen(moviefn,ii/pltfreq,[8,8],400);
    end
  end
end

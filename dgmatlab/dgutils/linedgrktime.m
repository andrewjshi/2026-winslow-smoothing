function u=linedgrktime(u,msh,data,phys,resname,dt,nsteps)
%LINEDGRKTIME  Line-DG Runge-Kutta Timestepping
%
%    Take nsteps with RK4 explicit timestepping scheme.
%
%    Syntax: u=linedgrktime(u,msh,data,phys,resname,dt,nsteps)
%
%    u:                Initial solution (input)
%    u:                Final solution (output)
%    msh, data, phys:  DG problem
%    resname:          Residual function
%    dt:               Stepsize
%    nsteps:           Total number of timesteps

if ~isfield(phys,'time')
  phys.time=0.0;
end
for ii=1:nsteps
  fprintf('.');

  if phys.viscous
      q = feval(resname, u, msh, data, phys, 'linedg_qeval');
      k1 = dt*feval(resname, u, q, msh, data, phys, 'linedg_assemble');
      phys.time = phys.time + dt/2;
      cu = addN(u, k1/2);
      q = feval(resname, cu, msh, data, phys, 'linedg_qeval');
      k2 = dt*feval(resname,cu, q, msh, data, phys, 'linedg_assemble');
      cu = addN(u, k2/2);
      q = feval(resname, cu, msh, data, phys, 'linedg_qeval');
      k3 = dt*feval(resname,cu, q, msh, data, phys, 'linedg_assemble');
      phys.time = phys.time + dt/2;
      cu = addN(u, k3);
      q = feval(resname, cu, msh, data, phys, 'linedg_qeval');
      k4 = dt*feval(resname, cu, q, msh, data, phys, 'linedg_assemble');
  else
      q0 = repmat(permute(0*u, [1,2,4,3]), [1,1,size(msh.p,1),1]);
      k1 = dt*feval(resname, u, q0, msh, data, phys, 'linedg_assemble');
      phys.time = phys.time + dt/2;
      cu = addN(u, k1/2);
      k2 = dt*feval(resname,cu, q0, msh, data, phys, 'linedg_assemble');
      cu = addN(u, k2/2);
      k3 = dt*feval(resname,cu, q0, msh, data, phys, 'linedg_assemble');
      phys.time = phys.time + dt/2;
      cu = addN(u, k3);
      k4 = dt*feval(resname,cu, q0, msh, data, phys, 'linedg_assemble');
  end
  u = addN(u, (k1 + 2*k2 + 2*k3 + k4) / 6);
  if any(isnan(u(:))), error('NaN in solution.'); end
end
fprintf('\n');

function c = addN(a, b)

c = a;
nb = size(b,2);
c(:,1:nb,:) = c(:,1:nb,:) + b;

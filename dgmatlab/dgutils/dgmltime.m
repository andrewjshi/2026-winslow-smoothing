function u=dgbdftime(u,msh,data,phys,jacname,odefun,tspan,tols,pltpars,moviefn)
%DGMLTIME  DG Timestepping using MATLAB's ODE solvers
%
%    Take nsteps with BDF implicit timestepping scheme.
%
%    NOTE: Currently only supports N==NU
%
%    Syntax: u=dgbdftime(u,msh,data,phys,jacname,odefun,tspan, ...
%                        tols,pltpars,moviefn)
%
%    u:                Initial solution (input)
%    u:                Final solution (output)
%    msh, data, phys:  DG problem
%    jacname:          Jacobian function
%    odefun:           ODE solver, 'ode15s' or 'ode45'
%    tspan:            Integration times and output times (see ML doc)
%    tols:             [abstol,reltol]
%
%    Optional plotting parameters:
%
%    pltpars:          Parameters to dgplot (cell array)
%    moviefn:          Directory and partial filename for images

if nargin<6 | isempty(odefun), odefun='ode15s'; end
if nargin<8 | isempty(tols),   tols=[1e-3,1e-2]; end
if nargin<9,  pltpars={}; end
if nargin<10, moviefn=''; end

M=mkmassmatrix(msh,data,size(u,2));
opts=odeset('abstol',tols(1),'reltol',tols(2),'outputfcn',@outfun, ...
            'mass',M,'mstatedependence','none', ...
            'masssingular','no','maxorder',3,'refine',1);
if strcmp(odefun,'ode15s')
  opts=odeset(opts,'jacobian',@jacfun);
end

pass.msh=msh;
pass.data=data;
pass.phys=phys;
pass.jacname=jacname;
pass.pltpars=pltpars;

[t,u]=feval(odefun,@resfun,tspan,u(:),opts,pass);
u=reshape(u(end,:),size(pass.msh.s,1),[],size(pass.msh.t,2));

function R=resfun(t,u,pass)

pass.phys.pars(1)=t;
u=reshape(u,size(pass.msh.s,1),[],size(pass.msh.t,2));
R=feval(pass.jacname,u,pass.msh,pass.data,pass.phys);
R=R(:);

function K=jacfun(t,u,pass)

pass.phys.pars(1)=t;
u=reshape(u,size(pass.msh.s,1),[],size(pass.msh.t,2));
[R,DK,OK]=feval(pass.jacname,u,pass.msh,pass.data,pass.phys);
K=mkjacobian(pass.msh.t2t,pass.data.egix,DK,OK);

function stop=outfun(t,u,cmd,pass)

switch cmd
 case 'done'
 case 'init'
 otherwise
  fprintf('Time = %g\n',t);
  u=reshape(u,size(pass.msh.s,1),[],size(pass.msh.t,2));
  dgplot(pass.msh,u,pass.pltpars{:});
end
stop=0;

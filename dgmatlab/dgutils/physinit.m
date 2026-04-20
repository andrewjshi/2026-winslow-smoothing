function phys=physinit(msh,bndcnds,flow,visc,vel0)

dim=size(msh.p,1);
if nargin<2, bndcnds=[]; end
phys.bndcnds=int32(bndcnds);

if nargin<3
  phys.uinf=[0];
  phys.pars=[];
  phys.viscous=true;
  return
end

if nargin<3 || isempty(flow), flow=[1.4,.1]; end
if nargin<4 || isempty(visc), visc=[100,0.72]; end
if nargin<5 || isempty(vel0), vel0=[1,zeros(1,dim-1)]; end

gamma=flow(1);
M0=flow(2);

phys.uinf=[1,vel0,1/gamma/M0^2/(gamma-1)+1/2];
phys.pars=visc;
phys.viscous=isnan(visc(1)) || ~isinf(visc(1));
phys.sensor_pars=zeros(1,4);
phys.time=0;

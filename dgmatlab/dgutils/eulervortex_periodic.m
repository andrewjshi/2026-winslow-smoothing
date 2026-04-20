function u=eulervortex_periodic(x,y,time,pars,geom,rep)

pars1=pars;
pars1(2)=0;
u0=eulervortex(x,y,time,pars1);
uinf=u0(1,:,1);

u=0*u0;
for i=-rep:rep
  for j=-rep:rep
    pars1=pars;
    pars1(5)=pars(5)-i*geom(1);
    pars1(6)=pars(6)-j*geom(2);
    u=u+eulervortex(x,y,time,pars1);
  end
end
u=u-((2*rep+1)^2-1)*repmat(uinf,[size(u,1),1,size(u,3)]);

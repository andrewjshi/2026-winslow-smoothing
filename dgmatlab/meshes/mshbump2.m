function msh=mshbump2(nref)
% BUMP2 (for Euler shock)

if nargin<1, nref=0; end

fd=inline('ddiff(drectangle(p,0,10,0,5),dcircle(p,3,-2.75,3))','p');
fdargs={};

phi0=asin(2.75/3.00);
s=block(8,-4)';
phi=phi0+(pi-phi0-phi0)*(1-s);
x=3+3*cos(phi);
y=-2.75+3*sin(phi);

pv=[0,0;x,y;10,0;10,5;0,5];
[p,t]=polymesh({pv},[1],[1,0],[4,1.6]);
for iref=1:nref
  [p,t]=uniref(p,t,1);
  p=bndproj(p,t,fd,fdargs{:});
end

bndexpr={'all(p(:,2)<1e-6)','all(p(:,1)>10-1e-6)', ...
         'all(p(:,2)>5-1e-6)','all(p(:,1)<1e-6)','true'};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,5);
msh=mshfix(msh);

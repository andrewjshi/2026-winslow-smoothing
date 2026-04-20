function msh=mshcollela(hmin,hmax,nref)

if nargin<1, hmin=0.05; end
if nargin<2, hmax=0.20; end
if nargin<3, nref=0; end

pv=[0,0; 0.6,0; 0.6,0.2-hmin; 0.6,0.2; 0.6+hmin,0.2;
    3,0.2; 3,1; 0,1];

[p,t]=polymesh({pv},[1],[1,0],[hmax,1.3]);
[p,t]=uniref(p,t,nref);

bndexpr={'all(p(:,1)<0+1e-6)','all(p(:,1)>3-1e-6)','true'};
msh=ml2msh(p,t,bndexpr,[],[]);
msh=mshcurved(msh,[]);
msh=mshfix(msh);

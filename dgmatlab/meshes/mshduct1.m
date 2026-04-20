function msh=mshduct1(mul1,quads)

if nargin<1, mul1=2; end
if nargin<2, quads=false; end

m=mul1+1;
n=3*mul1+1;

if quads
    [p,e,t]=qsquaremesh(n,m);
else
    [p,e,t]=squaremesh(n,m);
end

p(:,2)=p(:,2)*1;
p(:,1)=p(:,1)*3;

h=ductheight(p(:,1));
p(:,2)=h+(1-h).*p(:,2);

fd=inline('p(:,2)-ductheight(p(:,1))','p');
fdargs={};

bndexpr={'all(p(:,2)<1e-6) | all(p(:,1)>.1 & p(:,1)<2.9 & p(:,2)<1-1e-6)', ...
         'all(p(:,1)>3-1e-6)','all(p(:,2)>1-1e-6)','all(p(:,1)<1e-6)'};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,1);
msh=mshfix(msh);

function msh=mshcircincirc(m,n,hmin)

if nargin<1, m=21; end
if nargin<2, n=21; end
if nargin<3, hmin=0.1; end

rmin=1;
rmax=10;
phimin=0;
phimax=pi;

CT=block2tets(block([m,n],[1/hmin,1/hmin,1,1], ...
                    [rmin,phimin; rmax,phimin; rmin,phimax; rmax,phimax]));
XT=cat(1,CT(1,:,:).*cos(CT(2,:,:)),CT(1,:,:).*sin(CT(2,:,:)));
p=reshape(XT,2,[])';
t=reshape(1:numel(XT)/2,3,size(XT,3))';
[p,t]=fixmesh(p,t);

bndexpr={'all(p(:,2)<1e-6)', ...
         'all(sqrt(sum(p.^2,2))<1+1e-6)', ...
         'all(sqrt(sum(p.^2,2))>10-1e-6)'};
fd=inline('ddiff(dcircle(p,0,0,10),dcircle(p,0,0,1))','p');
fdargs={};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,[2,3]);
msh=mshfix(msh);

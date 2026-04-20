function msh=mshqcirchole(h1,h2,hparabola,width,height,growth)

if nargin<1, h1=0.4; end
if nargin<2, h2=1.5; end
if nargin<3, hparabola=[]; end % [-1.55,.2,.02,.5]
if nargin<4, width=4; end
if nargin<5, height=5; end
if nargin<6, growth=1.3; end

n1=ceil(pi/2/h1);
phi=linspace(pi,pi/2,n1+1)';
pv1=[cos(phi),sin(phi); 0,height; -width,height; -width,0];

[p,t]=polymesh({pv1},[1],[1,0],[h2,growth],@hh,hparabola);

fd=@dcircle;
fdargs={0,0,1};

bndexpr={'all(p(:,2)<1e-6)',...
         'all(sqrt(sum(p.^2,2))<1+1e-6)', ...
         'all(p(:,1)>-1e-6)', ...
         sprintf('all(p(:,2)>%g-1e-6)',height), ...
         sprintf('all(p(:,1)<-%g+1e-6)',width)};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,2);
msh=mshfix(msh);

function h=hh(p,hparabola)

h=inf+0*p(:,1);
if ~isempty(hparabola)
  p(:,1)=p(:,1)-hparabola(1);
  p=p(:,[2,1]);
  h1=hparabola(3)+hparabola(4)*abs(dparabola(p,hparabola(2)));
  h=min(h,h1);
end

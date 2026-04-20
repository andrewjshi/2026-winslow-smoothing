function msh=mshellipse1

fd=inline('ddiff(dcircle(p,0,0,5),dellipse(p,[.5,.25]))','p');
fdargs={};

n1=16;
n2=16;

phi1=2*pi*(1:n1)'/n1;
phi2=2*pi*(1:n2)'/n2;

x1=.5*cos(phi1);
x2=5*cos(phi2);
y1=0.25*sin(phi1);
y2=5*sin(phi2);

pv1=[x1,y1];
pv2=[x2,y2];

[p,t]=polymesh({pv1,pv2},[1,1],[0,1;1,0],[10,1.8]);
for ii=1:5
  p=bndproj(p,t,fd,fdargs{:});
end

bndexpr={'all((p(:,1)).^2+p(:,2).^2<2)','true'};
msh=ml2msh(p,t,bndexpr,fd,fdargs);
msh=mshcurved(msh,[1,2]);
msh=mshfix(msh);

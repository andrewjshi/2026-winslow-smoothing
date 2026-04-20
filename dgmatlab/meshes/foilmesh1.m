function [X,Y]=foilmesh1(geom,pars)

if nargin<1, geom='rae'; end
if nargin<2, pars=struct; end

% Default parameters
if ~isfield(pars,'NPNT'), pars.NPNT=[73,49,21]; end
if ~isfield(pars,'DS'), pars.DS=[1e-6,2e-5,1e-2]; end
if ~isfield(pars,'wderiv'), pars.wderiv=[5,5]; end
if ~isfield(pars,'angle'), pars.angle=0; end
if ~isfield(pars,'R'), pars.R=10; end
if ~isfield(pars,'wakeratio'), pars.wakeratio=200; end
if ~isfield(pars,'wakelength'), pars.wakelength=20; end
if ~isfield(pars,'extras'), pars.extras=[3,0,0]; end
if ~isfield(pars,'nsmooths'), pars.nsmooths=[1000,50]; end
if ~isfield(pars,'porder'), pars.porder=1; end

NF=pars.NPNT(1);
NY=pars.NPNT(2);
NW=pars.NPNT(3);
DS1=pars.DS(1);
DS2=pars.DS(2);
DS3=pars.DS(3);

switch geom
 case 'rae'
  load raebnd
  wakes=[41,217];
  if NF>0
    p=[x,y];
    wp=p(41:-1:1,:);
    p=p(41:217,:);
    phi=pars.angle*pi/180;
    p(:,1)=p(:,1)-1;
    p=p*[cos(phi),-sin(phi);sin(phi),cos(phi)];
    p(:,1)=p(:,1)+1;
    p(:,1)=p(:,1)/max(p(:,1));
    s=[0,cumsum(sqrt(sum(diff(p,[],1).^2,2)))'];
    s=s/max(s);
    ss=block(NF,-1/100,[0;1]); ss=ss/max(ss);
    
    for i=1:pars.extras(1)
      ss=[linspace(ss(1),ss(3),5),ss(4:end-3),linspace(ss(end-2),ss(end),5)];
    end
    
    for jj=1:pars.extras(2)
      [foo,ix]=min(abs(ss-0.775));
      ix=round((ix-1)/4)*4+1;
      ss=[ss(1:ix-5),linspace(ss(ix-4),ss(ix+4),17),ss(ix+5:end)];
    end
    ss=squant(ss,pars.porder);
   
    p=interp1(s,p,ss,'cubic');
    x=p(:,1); y=p(:,2);
    ws=[0,cumsum(sqrt(sum(diff(wp,[],1).^2,2)))'];
    ws=ws/max(ws);
    wss=block(NW,pars.wakeratio,[0;pars.wakelength/20])';
    
    for i=1:pars.extras(3)
      wss=[linspace(wss(1),wss(3),5)';wss(4:end)];
    end
    wss=squant(wss,pars.porder);
    
    wss=wss(2:end);
    w=interp1(ws,wp,wss,'cubic');
    p=[flipud(w); p; w];
    x=p(:,1); y=p(:,2);
    wakes=[NW,length(x)-NW+1];
  end
 case 'ht13'
  load ht13bnd
  p=[x,y];
  phi=pars.angle*pi/180;
  p(:,1)=p(:,1)-1;
  p=p*[cos(phi),-sin(phi);sin(phi),cos(phi)];
  p(:,1)=p(:,1)+1;
  s=[0,cumsum(sqrt(sum(diff(p,[],1).^2,2)))'];
  s=s/max(s);
  ss=block(NF,-1/100,[0;1]); ss=ss/max(ss);
  for i=1:4
    ss=[linspace(ss(1),ss(3),5),ss(4:end-3),linspace(ss(end-2),ss(end),5)];
  end
  ss=squant(ss,pars.porder);
  p=interp1(s,p,ss);
  x=p(:,1); y=p(:,2);
  w=block(NW,400,[1;pars.wakelength])';
  w=w(2:end);
  x=[flipud(w); x; w];
  y=[0*w; y; 0*w];
  wakes=[NW,length(x)-NW+1];
  x=flipud(x);
  y=flipud(y);
 otherwise
  error('Unknown foil geometry.');
end

p=[x,y];
s=[0;cumsum(sqrt(sum(diff(p,1,1).^2,2)))];
t=diff1d(p,s);
tlength=sqrt(sum(t.^2,2));
n=[-t(:,2)./tlength,t(:,1)./tlength];

p1=p;
n1=n;
NX=length(x);

for ismooth=1:pars.nsmooths(1)
  n1(2:end-1,:)=(n1(3:end,:)+n1(1:end-2,:))/2;
  n1length=sqrt(sum(n1.^2,2));
  n1=n1./[n1length,n1length];
end

p2=p1+pars.R*n1;

s2=[0;cumsum(sqrt(sum(diff(p2,1,1).^2,2)))];
t2=diff1d(p2,s2);
t2length=sqrt(sum(t2.^2,2));
n2=[-t2(:,2)./t2length,t2(:,1)./t2length];

n1=n;
for ismooth=1:pars.nsmooths(2)
  n1(2:end-1,:)=(n1(3:end,:)+n1(1:end-2,:))/2;
  n1length=sqrt(sum(n1.^2,2));
  n1=n1./[n1length,n1length];
end

ss=sdistr(NX,NY,pars.R,wakes,[DS1,DS2,DS3]);
ss=ss/pars.R;

% Linear TFI
%H0=(1-ss);
%H1=ss;
%X=p1(:,1*ones(1,NY)).*H0+p2(:,1*ones(1,NY)).*H1;
%Y=p1(:,2*ones(1,NY)).*H0+p2(:,2*ones(1,NY)).*H1;

% Hermite TFI
H0=2*ss.^3-3*ss.^2+1;
H1=ss.^3-2*ss.^2+ss;
H2=-2*ss.^3+3*ss.^2;
H3=ss.^3-ss.^2;
X=p1(:,1*ones(1,NY)).*H0+p2(:,1*ones(1,NY)).*H2+pars.wderiv(1)*n1(:,1*ones(1,NY)).*H1+pars.wderiv(2)*n2(:,1*ones(1,NY)).*H3;
Y=p1(:,2*ones(1,NY)).*H0+p2(:,2*ones(1,NY)).*H2+pars.wderiv(1)*n1(:,2*ones(1,NY)).*H1+pars.wderiv(2)*n2(:,2*ones(1,NY)).*H3;

if nargout==0
  % Plotting
  surf(X,Y,0*X),view(2),axis equal
  clear X Y
end

function s=squant(s,porder)

if porder==1, return; end
ns=length(s);
if mod(ns-1,porder)~=0
  error('Number of points %d not multiple of porder %d.',ns-1,porder);
end

for is=1:porder:ns-porder+1
  s(is:is+porder)=linspace(s(is),s(is+porder),porder+1);
end

function [X,Y]=foilmesh1(geom,pars)

% pars = [A,B,NY,DS1,DS2,DS3,angle,R]

if nargin<1, geom='rae'; end
if nargin<2, pars=[]; end
if length(pars)<1, NF=73; else NF=pars(1); end
if length(pars)<2, NY=49; else NY=pars(2); end
if length(pars)<3, NW=21; else NW=pars(3); end
if length(pars)<4, DS1=1e-6; else DS1=pars(4); end
if length(pars)<5, DS2=2e-5; else DS2=pars(5); end
if length(pars)<6, DS3=1e-2; else DS3=pars(6); end
if length(pars)<7, A=5; else A=pars(7); end
if length(pars)<8, B=5; else B=pars(8); end
if length(pars)<9, angle=0; else angle=pars(9); end
if length(pars)<10, R=10; else R=pars(10); end
if length(pars)<11, extras=0; else extras=pars(11); end
if length(pars)<13, nsmooths=[1000,50]; else nsmooths=pars(12:13); end
if length(pars)<14, extrasbnd=3; else extrasbnd=pars(14); end
if length(pars)<15, extrasmid=0; else extrasmid=pars(15); end
if length(pars)<16, extraswake=0; else extraswake=pars(16); end
if length(pars)<17, wakeratio=200; else wakeratio=pars(17); end
if length(pars)<18, wangle=0; else wangle=pars(18); end

switch geom
 case 'rae'
  load raebnd
  wakes=[41,217];
%  nsmooths=[1000,100];
  if NF>0
    p=[x,y];
    wp=p(41:-1:1,:);
    wp0=wp;
    wp(:,1)=wp(:,1)-wp0(1,1);
    wp(:,2)=wp(:,2)-wp0(1,2);
    phi=wangle*pi/180;
    wp=wp*[cos(phi),sin(phi);-sin(phi),cos(phi)];
    wp(:,1)=wp(:,1)+wp0(1,1);
    wp(:,2)=wp(:,2)+wp0(1,2);
    
    p=p(41:217,:);
    phi=angle*pi/180;
    p(:,1)=p(:,1)-1;
    p=p*[cos(phi),-sin(phi);sin(phi),cos(phi)];
    p(:,1)=p(:,1)+1;
    p(:,1)=p(:,1)/max(p(:,1));
    s=[0,cumsum(sqrt(sum(diff(p,[],1).^2,2)))'];
    s=s/max(s);
    ss=block(NF,-1/100,[0;1]); ss=ss/max(ss);

    for i=1:extrasbnd
      ss=[ss(1),mean(ss(1:2)),ss(2:end-1),mean(ss(end-1:end)),ss(end)];
    end

    ss=[ss(1:4),mean(ss(4:5)),ss(5:end)];
    
    for jj=1:extrasmid
      [foo,ix]=min(abs(ss-0.764));
      if ss(ix)>0.775
        ss=[ss(1:ix-1),(ss(ix-1)+ss(ix))/2,ss(ix:end)];
      else
        ss=[ss(1:ix),(ss(ix)+ss(ix+1))/2,ss(ix+1:end)];
      end
    end
  
    p=interp1(s,p,ss,'cubic');
    x=p(:,1); y=p(:,2);
    ws=[0,cumsum(sqrt(sum(diff(wp,[],1).^2,2)))'];
    ws=ws/max(ws);
    wss=block(NW,wakeratio,[0;1])';

    for i=1:extraswake
      wss=[wss(1);mean(wss(1:2));wss(2:end)];
    end
    
    wss=wss(2:end);
    w=interp1(ws,wp,wss,'cubic');
    p=[flipud(w); p; w];
    x=p(:,1); y=p(:,2);
    wakes=[NW,length(x)-NW+1];
  end
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

for ismooth=1:nsmooths(1)
  n1(2:end-1,:)=(n1(3:end,:)+n1(1:end-2,:))/2;
  n1length=sqrt(sum(n1.^2,2));
  n1=n1./[n1length,n1length];
end

p2=p1+R*n1;

s2=[0;cumsum(sqrt(sum(diff(p2,1,1).^2,2)))];
t2=diff1d(p2,s2);
t2length=sqrt(sum(t2.^2,2));
n2=[-t2(:,2)./t2length,t2(:,1)./t2length];

n1=n;
for ismooth=1:nsmooths(2)
  n1(2:end-1,:)=(n1(3:end,:)+n1(1:end-2,:))/2;
  n1length=sqrt(sum(n1.^2,2));
  n1=n1./[n1length,n1length];
end

ss=sdistr(NX,NY,R,wakes,[DS1,DS2,DS3]);
ss=ss/R;

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
X=p1(:,1*ones(1,NY)).*H0+p2(:,1*ones(1,NY)).*H2+A*n1(:,1*ones(1,NY)).*H1+B*n2(:,1*ones(1,NY)).*H3;
Y=p1(:,2*ones(1,NY)).*H0+p2(:,2*ones(1,NY)).*H2+A*n1(:,2*ones(1,NY)).*H1+B*n2(:,2*ones(1,NY)).*H3;

if nargout==0
  % Plotting
  surf(X,Y,0*X),view(2),axis equal
  clear X Y
end

function msh1=elliptic_smoothing_slide(msh,newp1,doplot,sliding_spec)
% ELLIPTIC_SMOOTHING_SLIDE  Variant of elliptic_smoothing that supports
% sliding the physical boundary along a smooth target curve. The original
% elliptic_smoothing already supports polygon-edge sliding via
% sliding_spec.corners; this variant adds:
%
%   sliding_spec.project_fn  - a function handle f(p)->q that projects
%                              a 1x2 point p onto the target boundary
%                              curve. Used after every PDE iteration to
%                              snap sliding nodes back onto the curve.
%   sliding_spec.pinned_xy   - Mx2 list of physical positions; boundary
%                              nodes within 1e-6 of any of these stay
%                              pinned (Dirichlet) for the duration.
%
% If sliding_spec is empty, behaves as fully pinned Winslow.

if nargin<2, newp1=msh.p1; msh=nodealloc(msh,msh.porder); end
if nargin<3, doplot=false; end
if nargin<4, sliding_spec=[]; end

D=size(msh.p,1);
msh1=msh;
msh.ncurved(:)=false;
msh.tcurved(:)=false;
msh.ecurved(:)=false;

data=dginit(msh);
e=getbndnodes(msh,data);
c=getcornernodes(msh,data);
ne=length(e);

[r,W,W1]=mapdg2cg(msh.p1);
N=size(W,1);

F=zeros(N,D);
newp=W*reshape(permute(newp1,[1,3,2]),[],size(newp1,2));

% Boundary classification.
slide_nodes   = [];
slide_proj_fn = [];
slide_mode    = 'none';
if ~isempty(sliding_spec) && isfield(sliding_spec,'project_fn')
    slide_mode    = 'curve';
    slide_proj_fn = sliding_spec.project_fn;
    pinned_xy = zeros(0,2);
    if isfield(sliding_spec,'pinned_xy')
        pinned_xy = sliding_spec.pinned_xy;
    end
    pin_tol = 1e-6;
    pinned_mask = false(length(e),1);
    for ii = 1:length(e)
        p = newp(e(ii),:);
        if ~isempty(pinned_xy) && min(vecnorm(pinned_xy - p, 2, 2)) < pin_tol
            pinned_mask(ii) = true;
        end
    end
    pinned_e = e(pinned_mask);
    bnd = [pinned_e, newp(pinned_e,:)];
    slide_nodes = e(~pinned_mask);
    fprintf('  Sliding Winslow (curve): %d pinned, %d sliding (of %d boundary nodes).\n', ...
        length(pinned_e), length(slide_nodes), length(e));
else
    bnd=[e,newp(e,:)];
end

u=W*reshape(permute(msh.p1,[1,3,2]),[],size(msh.p1,2));
for iter=1:100
  uold=u;
  u=sm1(msh,data,u,bnd);
  if strcmp(slide_mode,'curve')
    for si = 1:length(slide_nodes)
      ni = slide_nodes(si);
      u(ni,:) = slide_proj_fn(u(ni,:));
    end
  end
  uerror=norm(u(:)-uold(:));
  fprintf('  Smoothing, iteration %d: Norm = %16.12f\n',iter,uerror);
  if uerror<1e-7, break; end

  dgp=permute(reshape(W1'*u,size(msh.s,1),size(msh.t,2),[]),[1,3,2]);
  msh1.p1=dgp;

  if D==2 & doplot
    dgmeshplot_curved(msh1,2)
    line(u(:,1),u(:,2),'linesty','none','marker','.');
    drawnow
  end
end

if uerror>=1e-7, error('No convergence in smoothing.'); end

msh1.ncurved(:)=true;
msh1.tcurved(:)=true;
msh1.ecurved(:)=true;

J=geojac(msh);
J1=geojac(msh1);
fprintf('    Minimum Jacobian = %16.12f (initial) %16.12f (final)\n',min(J(:)),min(J1(:)));
if min(J1(:))<0, fprintf('  WARNING: Negative Jacobian\n'); end

function e=getbndnodes(msh,data)

ns=size(msh.s,1);
[j,it]=find(msh.t2t<0);
e=[];
for i=1:length(it)
  newe=ns*(it(i)-1)+double(data.egix(:,j(i)))+1;
  e=[e;newe];
end

evec=zeros(ns*size(msh.t,2),1);
evec(e)=1;
[r,foo,W]=mapdg2cg(msh.p1);
Wevec=W*evec;
e=find(Wevec>0);

function e=getcornernodes(msh,data)

D=size(msh.p,1);
ns=size(msh.s,1);
nt=size(msh.t,2);

ix=find(sum(msh.s<1e-3,2)==D);
e=ns*(repmat(1:nt,D+1,1)-1)+repmat(ix(:),1,nt);

evec=zeros(ns*size(msh.t,2),1);
evec(e)=1;
[r,foo,W]=mapdg2cg(msh.p1);
Wevec=W*evec;
e=find(Wevec>0);

function u=sm1(msh,data,u,bnd)

switch size(msh.p,1)
 case 2
  u=sm12d(msh,data,u,bnd);
 case 3
  u=sm13d(msh,data,u,bnd);
 otherwise
  error('Dimension not implemented.');
end

function u=sm12d(msh,data,u,bnd)

alpha=sm1alpha2d(msh,data,u,bnd);
alpha=cg2dg(msh,alpha);
u=cg2dg(msh,u);

p1x=permute(msh.p1(:,1,:),[1,3,2]);
p1y=permute(msh.p1(:,2,:),[1,3,2]);

phi=permute(data.gfs(:,1,:),[1,3,2]);
phiX=permute(data.gfs(:,2,:),[1,3,2]);
phiY=permute(data.gfs(:,3,:),[1,3,2]);

J11=phiX'*p1x;
J21=phiY'*p1x;
J12=phiX'*p1y;
J22=phiY'*p1y;
J=J11.*J22-J12.*J21;

iJ11=J22./J;
iJ22=J11./J;
iJ12=-J12./J;
iJ21=-J21./J;

ns=size(msh.s,1);
nt=size(msh.t,2);
KK=zeros(ns,ns,nt);
for it=1:nt
  phix=diag(iJ11(:,it))*phiX'+diag(iJ12(:,it))*phiY';
  phiy=diag(iJ21(:,it))*phiX'+diag(iJ22(:,it))*phiY';

  alpha1g=phi'*alpha(:,1,it);
  alpha2g=phi'*alpha(:,2,it);

  xxig=phix*u(:,1,it);
  xetag=phiy*u(:,1,it);
  yxig=phix*u(:,2,it);
  yetag=phiy*u(:,2,it);

  a=xetag.^2+yetag.^2;
  b=xxig.*xetag+yxig.*yetag;
  c=xxig.^2+yxig.^2;

  mul=data.gw.*J(:,it)/2;
  K=phix'*diag(mul.*a)*phix+phiy'*diag(mul.*c)*phiy- ...
    phix'*diag(mul.*b)*phiy-phiy'*diag(mul.*b)*phix;

  K=K+phi*diag(mul.*alpha1g)*phix+phi*diag(mul.*alpha2g)*phiy;

  KK(:,:,it)=K;
end

[ii,jj] = dgindices_mass(ns, nt);
[r,W]=mapdg2cg(msh.p1);
K=sparse(r(ii(:)),r(jj(:)),KK(:));

K(bnd(:,1),:)=0;
K(bnd(:,1),bnd(:,1))=speye(size(bnd,1),size(bnd,1));
F=zeros(size(K,1),2);
F(bnd(:,1),:)=bnd(:,[2,3]);
u=K\F;


function alpha=sm1alpha2d(msh,data,u,bnd)

u=cg2dg(msh,u);

p1x=permute(msh.p1(:,1,:),[1,3,2]);
p1y=permute(msh.p1(:,2,:),[1,3,2]);

phi=permute(data.gfs(:,1,:),[1,3,2]);
phiX=permute(data.gfs(:,2,:),[1,3,2]);
phiY=permute(data.gfs(:,3,:),[1,3,2]);

J11=phiX'*p1x;
J21=phiY'*p1x;
J12=phiX'*p1y;
J22=phiY'*p1y;
J=J11.*J22-J12.*J21;

iJ11=J22./J;
iJ22=J11./J;
iJ12=-J12./J;
iJ21=-J21./J;

ns=size(msh.s,1);
nt=size(msh.t,2);
RR1=zeros(ns,nt);
RR2=zeros(ns,nt);
for it=1:nt
  phix=diag(iJ11(:,it))*phiX'+diag(iJ12(:,it))*phiY';
  phiy=diag(iJ21(:,it))*phiX'+diag(iJ22(:,it))*phiY';

  xxig=phix*u(:,1,it);
  xetag=phiy*u(:,1,it);
  yxig=phix*u(:,2,it);
  yetag=phiy*u(:,2,it);

  E=xxig.^2+yxig.^2;
  F=xxig.*xetag+yxig.*yetag;
  G=xetag.^2+yetag.^2;

  g11=G;
  g22=E;
  g12=-F;
  g21=-F;

  mul=data.gw.*J(:,it)/2;
  R1=phix'*(mul.*g11)+phiy'*(mul.*g21);
  R2=phix'*(mul.*g12)+phiy'*(mul.*g22);
  RR1(:,it)=R1;
  RR2(:,it)=R2;
end

[r,foo,W]=mapdg2cg(msh.p1);
M = dgmass_ml(msh, data);
alpha=-(M\(W*[RR1(:),RR2(:)]));

function u=sm13d(msh,data,u,bnd)
error('elliptic_smoothing_slide: 3D not implemented.');

function u=cg2dg(msh,u)

[r,foo,W]=mapdg2cg(msh.p1);
u=permute(reshape(W'*u,size(msh.s,1),size(msh.t,2),[]),[1,3,2]);

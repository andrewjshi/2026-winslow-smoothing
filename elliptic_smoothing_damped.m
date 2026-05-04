function msh1=elliptic_smoothing_damped(msh,newp1,doplot,sliding_spec,opts)
% ELLIPTIC_SMOOTHING_DAMPED  Identical to elliptic_smoothing but adds a
% Picard damping factor and a configurable max-iter count.
%
%   opts.alpha    : damping factor in (0, 1]. u_{k+1} = (1-alpha)*u_k +
%                   alpha*sm1(u_k). alpha=1 reproduces elliptic_smoothing.
%                   Smaller alpha damps oscillations; typical 0.3-0.5.
%   opts.maxiter  : maximum Picard iterations (default 100).
%   opts.tol      : convergence tolerance on ||u_{k+1} - u_k|| (default 1e-7).

if nargin<2, newp1=msh.p1; msh=nodealloc(msh,msh.porder); end
if nargin<3, doplot=false; end
if nargin<4, sliding_spec=[]; end
if nargin<5, opts=struct(); end
if ~isfield(opts,'alpha'),   opts.alpha   = 0.5; end
if ~isfield(opts,'maxiter'), opts.maxiter = 100; end
if ~isfield(opts,'tol'),     opts.tol     = 1e-7; end

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

% Classify boundary nodes for sliding: if sliding_spec.corners is provided
% (domain corners in CCW order), boundary nodes near one of the corners
% stay pinned, all other boundary nodes are "sliding" --- they are freed
% from Dirichlet during each PDE solve and then projected back onto their
% assigned straight edge after the solve.
slide_nodes   = [];
slide_edge_ab = [];    % (M, 4) rows [ax ay bx by]
if ~isempty(sliding_spec) && isfield(sliding_spec,'corners')
    corners = sliding_spec.corners;
    K = size(corners,1);
    edges_ab = [corners, corners([2:K,1],:)];  % (K, 4)
    corner_tol = 1e-6;
    pinned_mask = false(length(e),1);
    slide_edge_idx = zeros(length(e),1);
    for ii = 1:length(e)
        p = newp(e(ii),:);
        if min(vecnorm(corners - p, 2, 2)) < corner_tol
            pinned_mask(ii) = true;
            continue;
        end
        best_ei = -1; best_d = inf;
        for ei = 1:K
            a = edges_ab(ei,1:2); b = edges_ab(ei,3:4);
            ab = b - a; ap = p - a;
            tt = max(0, min(1, (ap*ab')/(ab*ab')));
            d  = norm(p - (a + tt*ab));
            if d < best_d, best_d = d; best_ei = ei; end
        end
        slide_edge_idx(ii) = best_ei;
    end
    pinned_e = e(pinned_mask);
    bnd = [pinned_e, newp(pinned_e,:)];
    slide_nodes   = e(~pinned_mask);
    slide_edge_ab = edges_ab(slide_edge_idx(~pinned_mask),:);
    fprintf('  Sliding Winslow: %d pinned, %d sliding (of %d boundary nodes).\n', ...
        length(pinned_e), length(slide_nodes), length(e));
else
    bnd=[e,newp(e,:)];
end

u=W*reshape(permute(msh.p1,[1,3,2]),[],size(msh.p1,2));
fprintf('  Damped Picard: alpha=%g, maxiter=%d, tol=%.1e\n', ...
    opts.alpha, opts.maxiter, opts.tol);
for iter=1:opts.maxiter
  uold = u;
  u_new = sm1(msh,data,u,bnd);
  % Project sliding boundary nodes onto their assigned edges
  for si = 1:length(slide_nodes)
    ni = slide_nodes(si);
    a  = slide_edge_ab(si,1:2);
    b  = slide_edge_ab(si,3:4);
    ab = b - a; ap = u_new(ni,:) - a;
    tt = max(0, min(1, (ap*ab')/(ab*ab')));
    u_new(ni,:) = a + tt*ab;
  end
  u = (1 - opts.alpha) * uold + opts.alpha * u_new;
  uerror = norm(u(:) - uold(:));
  fprintf('  Smoothing, iteration %d: Norm = %16.12f\n',iter,uerror);
  if uerror<opts.tol, break; end

  dgp=permute(reshape(W1'*u,size(msh.s,1),size(msh.t,2),[]),[1,3,2]);
  msh1.p1=dgp;

  if D==2 & doplot
    dgmeshplot_curved(msh1,2)
    line(u(:,1),u(:,2),'linesty','none','marker','.');
    drawnow
  end
end

if uerror>=opts.tol
    fprintf('  WARNING: did not reach tol=%.1e (final residual %.4e). Returning last iterate.\n', ...
        opts.tol, uerror);
end
% Sync msh1.p1 to the final iterate (in case loop exited without the
% in-loop p1 update because the early-exit on tol skipped it).
msh1.p1 = permute(reshape(W1'*u, size(msh.s,1), size(msh.t,2), []), [1,3,2]);

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
% M=W*mkmassmatrix(msh,data,1)*W';
M = dgmass_ml(msh, data);
alpha=-(M\(W*[RR1(:),RR2(:)]));

function u=sm13d(msh,data,u,bnd)

alpha=sm1alpha3d(msh,data,u,bnd);
alpha=cg2dg(msh,alpha);
u=cg2dg(msh,u);

p1x=permute(msh.p1(:,1,:),[1,3,2]);
p1y=permute(msh.p1(:,2,:),[1,3,2]);
p1z=permute(msh.p1(:,3,:),[1,3,2]);

phi=permute(data.gfs(:,1,:),[1,3,2]);
phiX=permute(data.gfs(:,2,:),[1,3,2]);
phiY=permute(data.gfs(:,3,:),[1,3,2]);
phiZ=permute(data.gfs(:,4,:),[1,3,2]);

J11=phiX'*p1x;
J21=phiY'*p1x;
J31=phiZ'*p1x;
J12=phiX'*p1y;
J22=phiY'*p1y;
J32=phiZ'*p1y;
J13=phiX'*p1z;
J23=phiY'*p1z;
J33=phiZ'*p1z;
J=J11.*J22.*J33+J12.*J23.*J31+J13.*J21.*J32-J11.*J23.*J32-J12.*J21.*J33-J13.*J22.*J31;

iJ11=(J22.*J33-J23.*J32)./J;
iJ12=-(J12.*J33-J13.*J32)./J;
iJ13=(J12.*J23-J13.*J22)./J;
iJ21=-(J21.*J33-J23.*J31)./J;
iJ22=(J11.*J33-J13.*J31)./J;
iJ23=-(J11.*J23-J13.*J21)./J;
iJ31=(J21.*J32-J22.*J31)./J;
iJ32=-(J11.*J32-J12.*J31)./J;
iJ33=(J11.*J22-J12.*J21)./J;

ns=size(msh.s,1);
nt=size(msh.t,2);
KK=zeros(ns,ns,nt);
for it=1:nt
  phix=diag(iJ11(:,it))*phiX'+diag(iJ12(:,it))*phiY'+diag(iJ13(:,it))*phiZ';
  phiy=diag(iJ21(:,it))*phiX'+diag(iJ22(:,it))*phiY'+diag(iJ23(:,it))*phiZ';
  phiz=diag(iJ31(:,it))*phiX'+diag(iJ32(:,it))*phiY'+diag(iJ33(:,it))*phiZ';

  xxig=phix*u(:,1,it);
  xetag=phiy*u(:,1,it);
  xzetag=phiz*u(:,1,it);
  yxig=phix*u(:,2,it);
  yetag=phiy*u(:,2,it);
  yzetag=phiz*u(:,2,it);
  zxig=phix*u(:,3,it);
  zetag=phiy*u(:,3,it);
  zzetag=phiz*u(:,3,it);

  alpha1g=phi'*alpha(:,1,it);
  alpha2g=phi'*alpha(:,2,it);
  alpha3g=phi'*alpha(:,3,it);

  a11=xxig.^2+yxig.^2+zxig.^2;
  a22=xetag.^2+yetag.^2+zetag.^2;
  a33=xzetag.^2+yzetag.^2+zzetag.^2;
  a12=xxig.*xetag+yxig.*yetag+zxig.*zetag;
  a13=xxig.*xzetag+yxig.*yzetag+zxig.*zzetag;
  a23=xetag.*xzetag+yetag.*yzetag+zetag.*zzetag;

  b11=a22.*a33-a23.^2;
  b12=a13.*a23-a12.*a33;
  b13=a12.*a23-a13.*a22;
  b22=a11.*a33-a13.^2;
  b23=a13.*a12-a11.*a23;
  b33=a11.*a22-a12.^2;

  mul=data.gw.*J(:,it)/6;
  K=phix'*diag(mul.*b11)*phix+phix'*diag(mul.*b12)*phiy+phix'*diag(mul.*b13)*phiz+ ...
    phiy'*diag(mul.*b12)*phix+phiy'*diag(mul.*b22)*phiy+phiy'*diag(mul.*b23)*phiz+ ...
    phiz'*diag(mul.*b13)*phix+phiz'*diag(mul.*b23)*phiy+phiz'*diag(mul.*b33)*phiz;

  K=K+phi*diag(mul.*alpha1g)*phix+phi*diag(mul.*alpha2g)*phiy+phi*diag(mul.*alpha3g)*phiz;

  KK(:,:,it)=K;
end


[ii,jj] = dgindices_mass(ns, nt);
[r,W]=mapdg2cg(msh.p1);

rowix=r(ii(:));
colix=r(jj(:));
dat=KK(:);

bndix=ismember(rowix,bnd(:,1));
rowix=rowix(~bndix);
colix=colix(~bndix);
dat=dat(~bndix);

uniquebnd1=unique(bnd(:,1));
rowix=[rowix; uniquebnd1];
colix=[colix; uniquebnd1];
dat=[dat; ones(size(uniquebnd1))];

K=sparse(rowix,colix,dat);
F=zeros(size(K,1),3);
F(bnd(:,1),:)=bnd(:,[2,3,4]);
u=K\F;

function alpha=sm1alpha3d(msh,data,u,bnd)

u=cg2dg(msh,u);

p1x=permute(msh.p1(:,1,:),[1,3,2]);
p1y=permute(msh.p1(:,2,:),[1,3,2]);
p1z=permute(msh.p1(:,3,:),[1,3,2]);

phi=permute(data.gfs(:,1,:),[1,3,2]);
phiX=permute(data.gfs(:,2,:),[1,3,2]);
phiY=permute(data.gfs(:,3,:),[1,3,2]);
phiZ=permute(data.gfs(:,4,:),[1,3,2]);

J11=phiX'*p1x;
J21=phiY'*p1x;
J31=phiZ'*p1x;
J12=phiX'*p1y;
J22=phiY'*p1y;
J32=phiZ'*p1y;
J13=phiX'*p1z;
J23=phiY'*p1z;
J33=phiZ'*p1z;
J=J11.*J22.*J33+J12.*J23.*J31+J13.*J21.*J32-J11.*J23.*J32-J12.*J21.*J33-J13.*J22.*J31;

iJ11=(J22.*J33-J23.*J32)./J;
iJ12=-(J12.*J33-J13.*J32)./J;
iJ13=(J12.*J23-J13.*J22)./J;
iJ21=-(J21.*J33-J23.*J31)./J;
iJ22=(J11.*J33-J13.*J31)./J;
iJ23=-(J11.*J23-J13.*J21)./J;
iJ31=(J21.*J32-J22.*J31)./J;
iJ32=-(J11.*J32-J12.*J31)./J;
iJ33=(J11.*J22-J12.*J21)./J;

ns=size(msh.s,1);
nt=size(msh.t,2);
RR1=zeros(ns,nt);
RR2=zeros(ns,nt);
RR3=zeros(ns,nt);
for it=1:nt
  phix=diag(iJ11(:,it))*phiX'+diag(iJ12(:,it))*phiY'+diag(iJ13(:,it))*phiZ';
  phiy=diag(iJ21(:,it))*phiX'+diag(iJ22(:,it))*phiY'+diag(iJ23(:,it))*phiZ';
  phiz=diag(iJ31(:,it))*phiX'+diag(iJ32(:,it))*phiY'+diag(iJ33(:,it))*phiZ';

  xxig=phix*u(:,1,it);
  xetag=phiy*u(:,1,it);
  xzetag=phiz*u(:,1,it);
  yxig=phix*u(:,2,it);
  yetag=phiy*u(:,2,it);
  yzetag=phiz*u(:,2,it);
  zxig=phix*u(:,3,it);
  zetag=phiy*u(:,3,it);
  zzetag=phiz*u(:,3,it);

  a11=xxig.^2+yxig.^2+zxig.^2;
  a22=xetag.^2+yetag.^2+zetag.^2;
  a33=xzetag.^2+yzetag.^2+zzetag.^2;
  a12=xxig.*xetag+yxig.*yetag+zxig.*zetag;
  a13=xxig.*xzetag+yxig.*yzetag+zxig.*zzetag;
  a23=xetag.*xzetag+yetag.*yzetag+zetag.*zzetag;

  b11=a22.*a33-a23.^2;
  b12=a13.*a23-a12.*a33;
  b13=a12.*a23-a13.*a22;
  b22=a11.*a33-a13.^2;
  b23=a13.*a12-a11.*a23;
  b33=a11.*a22-a12.^2;

  mul=data.gw.*J(:,it)/6;
  R1=phix'*(mul.*b11)+phiy'*(mul.*b12)+phiz'*(mul.*b13);
  R2=phix'*(mul.*b12)+phiy'*(mul.*b22)+phiz'*(mul.*b23);
  R3=phix'*(mul.*b13)+phiy'*(mul.*b23)+phiz'*(mul.*b33);
  RR1(:,it)=R1;
  RR2(:,it)=R2;
  RR3(:,it)=R3;
end

[r,foo,W]=mapdg2cg(msh.p1);
%M=W*mkmassmatrix(msh,data,1)*W';
M = dgmass_ml(msh, data);
alpha=-(M\(W*[RR1(:),RR2(:),RR3(:)]));

function u=cg2dg(msh,u)

[r,foo,W]=mapdg2cg(msh.p1);
u=permute(reshape(W'*u,size(msh.s,1),size(msh.t,2),[]),[1,3,2]);

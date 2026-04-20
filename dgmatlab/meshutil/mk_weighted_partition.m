function part=mk_weighted_partition(npart,u,msh,data,phys,fname,dt,porder)

if nargin<6 | isempty(fname), fname=@dgnavierstokes; end
if nargin<7 | isempty(fname), dt=inf; end
if nargin<8, porder=msh.porder; end

t2t=msh.t2t'+1;

if porder<msh.porder
  u=dginterp(u,double(msh.porder),porder,size(msh.p,1));
  msh.tcurved(:)=false;
  msh=nodealloc(msh,porder);
  data=dginit(msh);
end

A=weight_matrix(u,msh,data,phys,fname,dt);
A=abs(A+A')/2;
ix=A~=0;
maxA=max(A(:));
mul=100/maxA;
A=ceil(A*mul);

nt=size(t2t,1);

part=metismex('PartGraphRecursive',A,npart,1)+1;
if all(part==2) % Seems like a bug in metis
  part(:)=1;
end

function C=weight_matrix(u,msh,data,phys,fname,dt)

if nargin<5, fname=@dgnavierstokes; end
if nargin<6, dt=inf; end

[rr,dd,oo]=feval(fname,u,msh,data,phys);
if ~isinf(dt)
  mm=mkmassmatrix(msh,data,size(rr,2),true);
  dd=mm-dt*dd;
  oo=-dt*oo;
end

ns=size(msh.s,1);
nc=size(rr,2);
nt=size(msh.t,2);
t2t=double(msh.t2t'+1);
t2n=double(bitand(uint32(msh.t2n),uint32(15))'+1);
t=double(msh.t'+1);
p=msh.p';
dim=size(p,2);
N=size(u,2);

egix = [];
for i=1:N
    egix = [egix; data.egix+(i-1)*ns+1];
end

C=zeros(nt,dim+1);
for i=1:nt
  ADinv=inv(dd(:,:,i));
  for jix=1:dim+1
    j=t2t(i,jix);
    if j>=1
      Aij=oo(:,:,jix,i);
      if j>i
        Aij=Aij';
        C(i,jix)=sqrt(sum(sum(abs(ADinv*Aij).^2)));
      else
        C(i,jix)=sqrt(sum(sum(abs(ADinv(:,egix(:,jix))*Aij).^2)));
      end
    elseif j>=1
      C(i,jix)=1e-3;
    end
  end
end

C=dualadjacency(t2t,C);

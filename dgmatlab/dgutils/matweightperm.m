function perm=matweightperm(u,msh,data,phys,fname,dt,GS,porder,partwise)

% Breadth-first

if nargin<5 | isempty(fname), fname=@dgnavierstokes; end
if nargin<6 | isempty(fname), dt=inf; end
if nargin<7 | isempty(fname), GS=false; end
if nargin<8, porder=msh.porder; end
if nargin<9, partwise=false; end

if porder>0
  error('This version only for porder=0');
end

if porder<msh.porder & msh.eltype==0
  u=dginterp(u,double(msh.porder),porder,size(msh.p,1),msh.eltype);
  msh.tcurved(:)=false;
  msh=nodealloc(msh,porder);
  data=dginit(msh);
end

if GS, error('GS not implemented in this version yet (easy to do).'); end

if partwise
  part=double(msh.tpart(1,:)'+1);
else
  part=ones(1,size(msh.t,2));
end

[rr,dd,oo]=feval(fname,u,msh,data,phys);
if ~isinf(dt)
  mm=mkmassmatrix(msh,data,size(rr,2),true);
  dd=mm-dt*dd;
  oo=-dt*oo;
end

ns=size(msh.s,1);
nc=size(rr,2);
nes = size(oo,1)/nc;
nt=size(msh.t,2);
nf=size(msh.t2t,1);
t2t=double(msh.t2t'+1);
t2n=double(bitand(uint32(msh.t2n),uint32(15))'+1);
t=double(msh.t'+1);
p=msh.p';
dim=size(p,2);

if msh.eltype==1
    if porder==0
        dd = reshape(mean(mean(reshape(dd, [ns,nc,ns,nc,nt]), ...
                               1), 3), [nc,nc,nt]);
        oo = reshape(mean(mean(reshape(oo, [nes,nc,ns,nc,nf,nt]), ...
                               1), 3), [nc,nc,nf,nt]);
    else
        error('foo');
    end
end

C=zeros(nt,nf);
for i=1:nt
  ADinv=inv(dd(:,:,i));
  for jix=1:nf
    j=t2t(i,jix);
    if j>=1 & part(i)==part(j)
      Aij=oo(:,:,jix,i);
      if j>i
        Aij=Aij';
      end
      C(i,jix)=sqrt(sum(sum(abs(ADinv*Aij).^2)));
    elseif j>=1 & part(i)~=part(j)
      C(i,jix)=1e-3;
    end
  end
end

w=zeros(nt,1);
for k=1:nt
  DC=zeros(nf,nf);
  for iix=1:nf
    i=t2t(k,iix);
    kix=t2n(k,iix);
    for jix=1:nf
      j=t2t(k,jix);
      if i>=1 & j>=1 & i~=j
        DC(iix,jix)=C(i,kix)*C(k,jix);
      end
    end
  end
  w(k)=sum(abs(DC(:)).^2);
end

perm=zeros(nt,1);
q=[];
for inbr=1:nt
  if isempty(q)
    [foo,i]=min(w);
    q=i;
  end
    
  [foo,ix]=min(w(q));
  i=q(ix); q(ix)=[];
  w=updiluw(w,t2t,t2n,C,i);
  perm(inbr)=i;

  for inix=1:nf
    in=t2t(i,inix);
    if in>0 & ~isinf(w(in)) & ~ismember(in,q)
      q=[q,in];
    end
  end
end

function w=updiluw(w,t2t,t2n,C,pi)

[nt,nf]=size(t2t);
w(pi)=inf;

for kix=1:nf
  k=t2t(pi,kix);
  if k>=1 & ~isinf(w(k))
    DC=zeros(nf,nf);
    for iiix=1:nf
      ii=t2t(k,iiix);
      kix=t2n(k,iiix);
      if ii>=1 & ~isinf(w(ii))
        for jjix=1:nf
          jj=t2t(k,jjix);
          if jj>=1 & ~isinf(w(jj))
            if ii~=jj
              DC(iiix,jjix)=C(ii,kix)*C(k,jjix);
            end
          end
        end
      end
    end
    w(k)=sum(abs(DC(:)).^2);
  end
end



% Sparse matrix version (slow)

% function perm=matweightperm(u,msh,data,phys,fname,dt,GS,porder,partwise)

% % Breadth-first

% if nargin<5 | isempty(fname), fname=@dgnavierstokes; end
% if nargin<6 | isempty(fname), dt=inf; end
% if nargin<7 | isempty(fname), GS=false; end
% if nargin<8, porder=msh.porder; end
% if nargin<9, partwise=false; end

% if porder<msh.porder
%   u=dginterp(u,double(msh.porder),porder,size(msh.p,1));
%   msh.tcurved(:)=false;
%   msh=nodealloc(msh,porder);
%   data=dginit(msh);
% end

% if GS, error('GS not implemented in this version yet (easy to do).'); end

% if partwise
%   part=double(msh.tpart(1,:)'+1);
% else
%   part=ones(1,size(msh.t,2));
% end

% [rr,dd,oo]=feval(fname,u,msh,data,phys);
% if ~isinf(dt)
%   mm=mkmassmatrix(msh,data,size(rr,2),true);
%   dd=mm-dt*dd;
%   oo=-dt*oo;
% end
% A=mkjacobian(msh.t2t,data.egix,dd,oo);

% ns=size(msh.s,1);
% nc=size(rr,2);
% nt=size(msh.t,2);
% t2t=double(msh.t2t'+1);
% t=double(msh.t'+1);
% p=msh.p';

% C=sparse(nt,nt);
% for i=1:nt
%   ix=ns*nc*(i-1)+(1:ns*nc);
%   ADinv=inv(full(A(ix,ix)));
%   for j=t2t(i,:)
%     if j>=1 & part(i)==part(j)
%       jx=ns*nc*(j-1)+(1:ns*nc);
%       C(i,j)=sqrt(sum(sum(abs(ADinv*A(ix,jx)).^2)));
%     elseif j>=1 & part(i)~=part(j)
%       C(i,j)=1e-3;
%     end
%   end
% end

% w=zeros(nt,1);
% for k=1:nt
%   DC=sparse(nt,nt);
%   for i=t2t(k,:)
%     for j=t2t(k,:)
%       if i>=1 & j>=1 & i~=j
%         DC(i,j)=C(i,k)*C(k,j);
%       end
%     end
%   end
%   w(k)=sum(abs(DC(:)).^2);
% end

% A=dualadjacency(t2t);
% perm=zeros(nt,1);

% q=[];
% for inbr=1:nt
%   if isempty(q)
%     [foo,i]=min(w);
%     q=i;
%   end
    
%   [foo,ix]=min(w(q));
%   i=q(ix); q(ix)=[];
%   w=updiluw(w,t2t,C,i);
%   perm(inbr)=i;
    
%   for in=find(A(:,i)')
%     if ~isinf(w(in))
%       q=[q,in];
%     end
%   end
% end

% function w=updiluw(w,t2t,C,pi)

% w(pi)=inf;
% nt=size(t2t,1);

% for k=t2t(pi,:)
%   if k>=1 & ~isinf(w(k))
%     DC=sparse(nt,nt);
%     for ii=t2t(k,:)
%       if ii>=1 & ~isinf(w(ii))
%         for jj=t2t(k,:)
%           if jj>=1 & ~isinf(w(jj))
%             if ii~=jj
%               DC(ii,jj)=C(ii,k)*C(k,jj);
%             end
%           end
%         end
%       end
%     end
%     w(k)=sum(abs(DC(:)).^2);
%   end
% end




%%% Original version, just pick minimum

% function p=matweightperm_orig(u,msh,data,phys,fname,dt,GS,porder)

% if nargin<5 | isempty(fname), fname=@dgnavierstokes; end
% if nargin<6 | isempty(fname), dt=inf; end
% if nargin<7 | isempty(fname), GS=false; end
% if nargin<8, porder=msh.porder; end

% if porder<msh.porder
%   u=dginterp(u,double(msh.porder),porder,size(msh.p,1));
%   msh.tcurved(:)=false;
%   msh=nodealloc(msh,porder);
%   data=dginit(msh);
% end

% [rr,dd,oo]=feval(fname,u,msh,data,phys);
% if ~isinf(dt)
%   mm=mkmassmatrix(msh,data,size(rr,2),true);
%   dd=mm-dt*dd;
%   oo=-dt*oo;
% end
% A=mkjacobian(msh.t2t,data.egix,dd,oo);

% ns=size(msh.s,1);
% nc=size(rr,2);
% nt=size(msh.t,2);
% t2t=double(msh.t2t'+1);

% C=sparse(nt,nt);
% for i=1:nt
%   ix=ns*nc*(i-1)+(1:ns*nc);
%   ADinv=inv(full(A(ix,ix)));
%   for j=t2t(i,:)
%     if j>=1
%       jx=ns*nc*(j-1)+(1:ns*nc);
%       C(i,j)=sqrt(sum(sum(abs(ADinv*A(ix,jx)).^2)));
%     end
%   end
% end

% w=zeros(nt,1);
% for k=1:nt
%   DC=sparse(nt,nt);
%   if GS
%     for i=t2t(k,:)
%       if i>=1
%         DC(i,1)=C(i,k);
%       end
%     end
%   else
%     for i=t2t(k,:)
%       for j=t2t(k,:)
%         if i>=1 & j>=1 & i~=j
%           DC(i,j)=C(i,k)*C(k,j);
%         end
%       end
%     end
%   end
%   w(k)=sum(abs(DC(:)).^2);
% end

% p=zeros(1,nt);
% for i=1:nt
%   [foo,p(i)]=min(w);
%   w(p(i))=inf;
  
%   for k=t2t(p(i),:)
%     if k>=1 & ~isinf(w(k))
%       DC=sparse(nt,nt);
%       if GS
%         for ii=t2t(k,:)
%           if ii>=1 & ~isinf(w(ii))
%             DC(ii,1)=C(ii,k);
%           end
%         end
%       else
%         for ii=t2t(k,:)
%          if ii>=1 & ~isinf(w(ii))
%            for jj=t2t(k,:)
%              if jj>=1 & ~isinf(w(jj))
%                if ii~=jj
%                  DC(ii,jj)=C(ii,k)*C(k,jj);
%                end
%              end
%            end
%          end
%         end
%       end
%       w(k)=sum(abs(DC(:)).^2);
%     end
%   end
% end

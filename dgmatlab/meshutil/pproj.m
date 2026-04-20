function p=pproj(p,d0,fd,varargin)

deps=sqrt(eps)*max(max(p)-min(p));

if size(p,2)==2
  d=feval(fd,p,varargin{:});
  dgradx=(feval(fd,[p(:,1)+deps,p(:,2)],varargin{:})-d)/deps;
  dgrady=(feval(fd,[p(:,1),p(:,2)+deps],varargin{:})-d)/deps;
  dgrad2=dgradx.^2+dgrady.^2;
  dgrad2(dgrad2==0)=1;
  dd=d-d0;
  p=p-[dd.*dgradx./dgrad2,dd.*dgrady./dgrad2];
elseif size(p,2)==3
  d=feval(fd,p,varargin{:});
  dgradx=(feval(fd,[p(:,1)+deps,p(:,2),p(:,3)],varargin{:})-d)/deps;
  dgrady=(feval(fd,[p(:,1),p(:,2)+deps,p(:,3)],varargin{:})-d)/deps;
  dgradz=(feval(fd,[p(:,1),p(:,2),p(:,3)+deps],varargin{:})-d)/deps;
  dgrad2=dgradx.^2+dgrady.^2+dgradz.^2;
  dgrad2(dgrad2==0)=1;
  dd=d-d0;
  p=p-[dd.*dgradx./dgrad2,dd.*dgrady./dgrad2,dd.*dgradz./dgrad2];
end

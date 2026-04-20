function expr=dgeval(u,utype,dim)

if nargin<3, dim=size(u,2)-2; end
gamma=1.4;

if isnumeric(utype)
  expr=u(:,utype,:);
elseif char(utype)
  switch utype
   case 'r'
    expr=u(:,1,:);
   case 'ru'
    expr=u(:,2,:);
   case 'rv'
    expr=u(:,3,:);
   case 'rw'
    expr=u(:,4,:);
   case 'rE'
    expr=u(:,dim+2,:);
   case 'u'
    expr=u(:,2,:)./u(:,1,:);
   case 'v'
    expr=u(:,3,:)./u(:,1,:);
   case 'w'
    expr=u(:,4,:)./u(:,1,:);
   case 'E'
    expr=u(:,dim+2,:)./u(:,1,:);
   case 'ru2'
    expr=sum(u(:,2:dim+1,:).^2,2);
   case 'ruabs'
    expr=sqrt(dgeval(u,'ru2',dim));
   case 'uabs'
    expr=dgeval(u,'ruabs',dim)./u(:,1,:);
   case 'p'
    expr=(1.4-1)*(u(:,dim+2,:)-1/2*(dgeval(u,'ru2',dim))./u(:,1,:));
   case 's'
    expr=dgeval(u,'p',dim)./u(:,1,:).^1.4;
   case 'T'
    expr=dgeval(u,'p',dim)./u(:,1,:)/(gamma-1);
   case 'M'
    expr=sqrt(dgeval(u,'ru2',dim)./u(:,1,:)./dgeval(u,'p',dim)/1.4);
   otherwise
    error('Unknown expression.');
  end
end

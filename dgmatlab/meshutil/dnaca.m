function d=dnaca(p,a,theta,varargin)

p=protate(p,-theta);
e1=polyval([a(5:-1:2),0],p(:,1));
e2=a(1)^2*p(:,1);
d1=(p(:,2)-e1).^2-e2;
d2=(p(:,2)+e1).^2-e2;
d=max(d1,d2);

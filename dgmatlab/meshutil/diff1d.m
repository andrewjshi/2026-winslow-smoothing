function dy=diff1d(y,s)

s=s(:);
dy=0*y;
dy(2:end-1,:)=(y(3:end,:)-y(1:end-2,:))./repmat(s(3:end)-s(1:end-2),1,size(y,2));
dy(1,:)=(y(2,:)-y(1,:))/(s(2)-s(1));
dy(end,:)=(y(end,:)-y(end-1,:))/(s(end)-s(end-1));

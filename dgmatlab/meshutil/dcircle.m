function d=dcirc(p,xc,yc,r)

d=sqrt((p(:,1)-xc).^2+(p(:,2)-yc).^2)-r;

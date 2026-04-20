function d=drect(p,x1,x2,y1,y2)

d=-min(min(min(-y1+p(:,2),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));

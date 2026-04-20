function J=geojac(msh,data)

if nargin<2, data=dginit(msh); end

dim=size(msh.p,1);
if dim==1
  p1x=permute(msh.p1(:,1,:),[1,3,2]);
  phiX=permute(data.gfs(:,2,:),[1,3,2]);
  
  J=phiX'*p1x;
elseif dim==2
  p1x=permute(msh.p1(:,1,:),[1,3,2]);
  p1y=permute(msh.p1(:,2,:),[1,3,2]);
  
  phiX=permute(data.gfs(:,2,:),[1,3,2]);
  phiY=permute(data.gfs(:,3,:),[1,3,2]);
  
  xX=phiX'*p1x;
  xY=phiY'*p1x;
  yX=phiX'*p1y;
  yY=phiY'*p1y;
  
  J=xX.*yY-xY.*yX;
elseif dim==3
  p1x=permute(msh.p1(:,1,:),[1,3,2]);
  p1y=permute(msh.p1(:,2,:),[1,3,2]);
  p1z=permute(msh.p1(:,3,:),[1,3,2]);
  
  phiX=permute(data.gfs(:,2,:),[1,3,2]);
  phiY=permute(data.gfs(:,3,:),[1,3,2]);
  phiZ=permute(data.gfs(:,4,:),[1,3,2]);
  
  xX=phiX'*p1x;
  xY=phiY'*p1x;
  xZ=phiZ'*p1x;
  yX=phiX'*p1y;
  yY=phiY'*p1y;
  yZ=phiZ'*p1y;
  zX=phiX'*p1z;
  zY=phiY'*p1z;
  zZ=phiZ'*p1z;
  
  J=xX.*yY.*zZ+xY.*yZ.*zX+xZ.*yX.*zY-xX.*yZ.*zY-xY.*yX.*zZ-xZ.*yY.*zX;
end

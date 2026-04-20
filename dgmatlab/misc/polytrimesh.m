function [p,t] = polytrimesh(pvs,holes,szfcn,cmd)
%POLYTRIMESH Unstructured meshing of polygonal regions
%
%   Syntax: [P,T] = POLYTRIMESH(PVS,HOLES,CMD,SZFCN)
%
%   PVS: Cell array of polygons (Nx2 arrays)
%   HOLES: x,y coordinates of holes in geometry (Mx2 array)
%   SZFCN: Size function, expression (string, default 'inf')
%   CMD: Command to triangle (string, default 'puq28.6Q')
%
%   Example:
%      pv1 = [0,0;1,0;1,1;.5,.5;0,1;0,0];
%      pv2 = [.2,.2;.4,.2;.4,.4;.39,.4;.2,.4;.2,.2];
%      holes = [.3,.3];
%      szfcn = '0.01+0.3*abs(sqrt((x-.5)^2+(y-.5)^2)-.5)';
%      [p,t]=polytrimesh({pv1,pv2}, holes, szfcn);
%      simpplot(p,t)

if ~iscell(pvs), pvs={pvs}; end
if nargin<3 | isempty(szfcn), szfcn='inf'; end
if nargin<4 | isempty(cmd), cmd='puq28.6Q'; end

pv=[];
seg=[];
for i=1:numel(pvs)
  cpv=pvs{i};
  closed=false;
  if size(cpv,1)>1 & norm(cpv(1,:)-cpv(end,:))<1e-10
    closed=true;
    cpv=cpv(1:end-1,:);
  end
  
  if size(cpv,1)>1
    if closed
      cseg=[1:size(cpv,1); 2:size(cpv,1),1];
    else
      cseg=[1:size(cpv,1)-1; 2:size(cpv,1)];
    end
    
    seg=[seg, cseg+size(pv,1)];
  end
  
  pv=[pv;cpv];
end

[p,x,x,t,x,x,x,x,x]=trianglemex(szfcn, cmd, pv',[],[],[],[],[], ...
                                int32(seg),[],holes',[]);
p=p';
t=t';

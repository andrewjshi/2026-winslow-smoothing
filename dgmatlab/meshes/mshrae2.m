function [msh]=mshrae2(porder,resolution)

if porder>=3
    error('maximum value of porder is 2');
end

load rae_mesh;
nn=round((33-1)/resolution)+1; 
nq=round((nn-1)/porder)+1;
x=rae(1:resolution:end,1:resolution:end,1);
y=rae(1:resolution:end,1:resolution:end,2);
[nx,ny]=size(x);
xp=x(1:porder:nx,1:porder:ny);
yp=y(1:porder:nx,1:porder:ny);
[na,nr]=size(xp);

% Create p
p=[xp(:) yp(:)];
p(1:nq,:)=[]; % Remove the duplicate nodes

% Create t
t0=([1,2,na+2; 1,na+2,na+1]); 
t=[];
for i=0:nr-2
    for j=0:na-2
        t=[t;t0+j+i*(na)];
    end
end

% Offset the elements for duplicate nodes
t=t-nq;
for i=1:nq-1
    g=2*i-1;
    t(g,1)=na-nq-i+1;
    t(g,2)=na-nq-i;
    t(g+1,1)=na-nq-i+1;
end
t(2*nq-1,1)=na-2*nq+1;
t(2*nq-1,2)=1;
t(2*nq,1)=na-2*nq+1;


% Create msh
bndexpr={'all(sqrt(sum(p.^2,2))<2)','true'};
msh=ml2msh(p,t,bndexpr);
msh.ncurved=logical(ones(1,size(msh.p,2)));
msh.ecurved=logical(ones(3,size(msh.t,2)));
msh.tcurved=logical(ones(1,size(msh.t,2)));
[msh.s,msh.tlocal]=lagrangepnts(porder,2);
[msh.sbnd,msh.tbndlocal]=lagrangepnts(porder,1);

ns=size(msh.s,1);
nt=size(msh.t,2);
msh.porder=int32(porder);
dim=2;

% Create msh.p1
if porder==1
    msh.p1=zeros(ns,dim,nt);
    for comp=1:dim
      for trinode=1:dim+1
        dp=msh.s(:,trinode)*msh.p(comp,double(msh.t(trinode,:))+1);
        msh.p1(:,comp,:)=msh.p1(:,comp,:)+permute(dp,[1,3,2]);
      end
    end
else
    p1=[x(:) y(:)];     
    msh.p1=zeros(ns,dim,nt); k=1;
    for i=1:porder:ny-porder
        for j=1:porder:nx-porder
            g = (i-1)*nx+j;
            ind1=index1(g,porder,nx);
            msh.p1(:,:,k)=p1(ind1,:);
            ind2=index2(g,porder,nx);
            msh.p1(:,:,k+1)=p1(ind2,:);
            k=k+2;
        end
    end
end

p1=[x(:) y(:)]; 
ps=p1(nn:(369-1)/resolution-nn+1,:);
% figure(1);plot(ps(:,1),ps(:,2),'*');
% pause
np=size(ps,1);
msh.dd=zeros(ns,1,nt);
for ii=1:nt    
    for jj=1:ns        
        ts=repmat(msh.p1(jj,:,ii),[np 1]);
        ds=sqrt(sum((ts-ps).^2,2));                
        msh.dd(jj,1,ii)=min(ds);
    end
end
    
function ind1=index1(g,porder,nx)
if porder==2
    ind1=g+[2*nx+2 nx+1 0 nx+2 1 2];    
end

function ind2=index2(g,porder,nx)
if porder==2
    ind2=g+[2*nx nx 0 2*nx+1 nx+1 2*nx+2];    
end



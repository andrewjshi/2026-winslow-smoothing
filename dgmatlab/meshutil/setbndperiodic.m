function [t2t,t2n]=setbndperiodic(p0,t,t2t,t2n,map)
%SETBNDPERIODIC  Set periodic boundaries in t2t,t2n matrices.
%
%    Syntax: [t2t,t2n]=setbndperiodic(p0,t,t2t,t2n,map)
%
%    p0,t,t2t,t2n:
%          mesh in standard MATLAB format, boundaries 1,2,3,...
%          correspond to the numbers 0,-1,-2,... in t2t
%
%    map:  describes which boundaries are periodic:
%              map={bnd1,expr1,bnd2,expr2;
%                   ....,.....,....,.....}
%          will make bnd1 and bnd2 periodic, where elements
%          are matched based on the expressions in expr1
%          (expressions depending on nodes p)
%
%          Ex: map={1,'p(:,1)',2,'2*p(:,2)'} will make bnds 1,2
%          periodic, where x on bnd 1 is matched to 2*x on bnd 2

dim=size(p0,2);
nf=size(t2t,2);
%if ~dimnv2eltype(dim,size(t,2))==t_simplex
%    error('only implemented for simplex elements');
%end

for ii=1:size(map,1)
  [val1,elem1,face1,six1]=getval(map{ii,1},map{ii,2},p0,t,t2t);
  [val2,elem2,face2,six2]=getval(map{ii,3},map{ii,4},p0,t,t2t);
  nelem=length(elem1);

  [foo,foo,ix]=unique([val1;val2],'rows');
  ix=reshape(ix,nelem,2);
  for ii=1:nelem
    jj=find(ix(ii,1)==ix(:,2));
    t2t(elem1(ii),face1(ii))=elem2(jj);
    t2n(elem1(ii),face1(ii))=face2(jj)+2^7+2^4*(permmatch(six1(ii,:),six2(jj,:),dim)-1);
    t2t(elem2(jj),face2(jj))=elem1(ii);
    t2n(elem2(jj),face2(jj))=face1(ii)+2^7+2^4*(permmatch(six2(jj,:),six1(ii,:),dim)-1);
  end
end

function p=permmatch(n1,n2,dim)

switch dim
 case 1
  p=1;
 case 2
  p=1;
 case 3
   p=n2(find(n1==1));
   %p=find(n2==n1(1));
end

function x=snap(x)

tol=sqrt(eps);
x=tol*round(x/tol);

function [val,elem,face,six]=getval(bnd,expr,p0,t,t2t)

eltype=dimnv2eltype(size(p0,2),size(t,2));
val=[];
six=[];
[elem,face]=find(t2t==-bnd+1);
map=mkfacemap(size(p0,2),eltype);
for ii=1:length(elem)
  tix=map(:,face(ii))';
  n=t(elem(ii),tix);
  
  p=p0(n,:);
  v=eval(expr);
  
  [val(:,:,ii),six(:,ii)]=sortrows(snap(v));
end
val=reshape(permute(val,[3,1,2]),ii,[]);
six=six';

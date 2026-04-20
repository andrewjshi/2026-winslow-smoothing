function msh1=mshtri2prism(msh,nlayer,dy,comp,bndexpr,periodic)

icomp=setdiff(1:3,comp);
if periodic
  str=sprintf('p(:,[%d,%d])',icomp(1),icomp(2));
  periodic_map={1,str,2,str};
else
  periodic_map={};
end

p=msh.p';
t=msh.t'+1;
np=size(p,1);
nt=size(t,1);
if mod(comp,2)==0
  t=t(:,[2,1,3]);
end

pp=zeros([size(p,1),3,nlayer+1]);
for i=1:nlayer+1
  pp(:,icomp,i)=p;
  pp(:,comp,i)=(i-1)*dy;
end
pp=reshape(permute(pp,[1,3,2]),[],3);

tt=zeros([size(t,1),4,3,nlayer]);
for i=1:nlayer
  for it=1:nt
    t1=t(it,:)+(i-1)*np;
    t2=t(it,:)+i*np;
    tt(it,:,:,i)=gen_split_prism_layer(t1,t2)';
  end
end
tt=reshape(permute(tt,[3,1,4,2]),[],4);

bndexpr={sprintf('p(:,%d)<1e-6',comp), ...
         sprintf('p(:,%d)>%f-1e-6',comp,nlayer*dy), ...
         bndexpr{:}};

msh1=ml2msh(pp,tt,bndexpr,[],[],periodic_map);
msh1=mshbndreorder(msh1,[3:length(bndexpr),1,2]);
msh1=mshcurved(msh1,[]);
msh1.ncurved(:)=true; msh1.tcurved(:)=true; msh1.ecurved(:)=true;
msh1=nodealloc(msh1,msh.porder);

msh0=msh;
msh0.tcurved(:)=false; msh0.ncurved(:)=false; msh0.ecurved(:)=false;
msh0=nodealloc(msh0,msh.porder);

itet=0;
for i=1:nlayer
  for it=1:nt
    for j=1:3
      itet=itet+1;
      p1=msh0.p1(:,:,it);
      p2=msh.p1(:,:,it);
      p3=msh1.p1(:,:,itet);
      p4=p3;
      p1=round(1e8*p1)/1e8;
      p3=round(1e8*p3)/1e8;
      [foo,ix]=ismember(p3(:,icomp),p1,'rows');
      p4(:,icomp)=p2(ix,:);
      msh1.p1(:,:,itet)=p4;
    end
  end
end

function ts=gen_split_prism_layer(t1,t2)

pr=[t1,t2];
[foo,ix]=min(pr);
map1=[1,2,3,4,5,6;
      2,3,1,5,6,4;
      3,1,2,6,4,5];
pr=pr(map1(ix,:));
if pr(2)<pr(3)
  ts=[1,2,3,6;
      1,2,6,5;
      1,5,6,4];
else
  ts=[1,2,3,5;
      1,5,3,6;
      1,5,6,4];
end
ts=pr(ts);

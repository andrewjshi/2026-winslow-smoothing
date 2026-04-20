function [msh1,u1]=mshuquad2hex(msh,u,nlayer,dy,comp,bndexpr,periodic)

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
    t=t(:,[2,1,4,3]);
end

pp=zeros([size(p,1),3,nlayer+1]);
for i=1:nlayer+1
  pp(:,icomp,i)=p;
  pp(:,comp,i)=(i-1)*dy;
end
pp=reshape(permute(pp,[1,3,2]),[],3);

tt=zeros([size(t,1),8,nlayer]);
for i=1:nlayer
  for it=1:nt
    t1=t(it,:)+(i-1)*np;
    t2=t(it,:)+i*np;
    tt(it,:,i)=gen_hex_layer(t1,t2)';
  end
end
tt=reshape(permute(tt,[1,3,2]),[],8);

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

u1 = zeros(size(msh1.p1,1), size(u,2), size(msh1.p1,3));

ihex=0;
for i=1:nlayer
    for it=1:nt
        ihex=ihex+1;
        p1=msh0.p1(:,:,it);
        p2=msh.p1(:,:,it);
        p3=msh1.p1(:,:,ihex);
        p4=p3;
        p1=round(1e8*p1)/1e8;
        p3=round(1e8*p3)/1e8;
        [foo,ix]=ismember(p3(:,icomp),p1,'rows');
        p4(:,icomp)=p2(ix,:);
        u1(:,:,ihex) = u(ix,:,it);
        msh1.p1(:,:,ihex)=p4;
    end
end

function ts=gen_hex_layer(t1,t2)

ts=[t1,t2];

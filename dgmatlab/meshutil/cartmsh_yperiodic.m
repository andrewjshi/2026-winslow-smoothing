function msh=cartmsh_yperiodic(msh)

y=msh.p1(:,2,:);
ymin=min(y(:));
ymax=max(y(:));

bndexpr=msh.bndexpr;
bndexpr={sprintf('p(:,2)<%.12g',ymin+1e-6),sprintf('p(:,2)>%.12g',ymax-1e-6),bndexpr{:}};
periodic_map={1,'p(:,[1,3])',2,'p(:,[1,3])'};

p=msh.p';
t=double(msh.t'+1);
t2t=double(msh.t2t'+1);
t2n=double(msh.t2n'+1);

t2t(t2t<0)=0;
t2t=setbndnbrs(p,t,t2t,bndexpr);
[t2t,t2n]=setbndperiodic(p,t,t2t,t2n,periodic_map);

msh.t2t=int32(t2t'-1);
msh.t2n=int32(t2n'-1);
msh.t2t(msh.t2t<0)=msh.t2t(msh.t2t<0)+2;

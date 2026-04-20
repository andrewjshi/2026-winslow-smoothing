function fro=fro2ml(fname)

fro=struct;
fid=fopen(fname,'r');

ns=fscanf(fid,'%d',6);
ne=ns(1); np=ns(2);
ncv=ns(5); nsf=ns(6);
fgetl(fid);

p=zeros(np,3);
for ii=1:np
  jp=fscanf(fid,'%d',1);
  p(ii,:)=fscanf(fid,'%f',3)';
  fgetl(fid);
end
fro.p=p;

e=zeros(ne,4);
for ii=1:ne
  je=fscanf(fid,'%d',1);
  e(ii,:)=fscanf(fid,'%d',4)';
  fgetl(fid);
end
fro.e=e;

pc={};
for ii=1:ncv
  npc=fscanf(fid,'%d',2); npc=npc(2);
  pc{ii}=zeros(npc,2);
  for jj=1:npc
    pc{ii}(jj,:)=fscanf(fid,'%f',2)';
  end
  fgetl(fid);
end
fro.pc=pc;

ps={};
for ii=1:nsf
  nps=fscanf(fid,'%d',2); nps=nps(2);
  ps{ii}=zeros(nps,3);
  for jj=1:nps
    ps{ii}(jj,:)=fscanf(fid,'%f',3)';
  end
  fgetl(fid);
end
fro.ps=ps;

u={};
for ii=1:ncv
  nu=fscanf(fid,'%d',2); nu=nu(2);
  u{ii}=zeros(2,3,nu);
  for jj=1:nu
    u{ii}(:,:,jj)=reshape(fscanf(fid,'%f',6),[3,2])';
  end
  fgetl(fid);
end
fro.u=u;

vw={};
for ii=1:nsf
  nvw=fscanf(fid,'%d',3); nv=nvw(2); nw=nvw(3);
  vw{ii}=zeros(4,3,nv,nw);
  for jj=1:nw
    for kk=1:nv
      vw{ii}(:,:,kk,jj)=reshape(fscanf(fid,'%f',12),[3,4])';
    end
  end
  fgetl(fid);
end
fro.vw=vw;

fclose(fid);

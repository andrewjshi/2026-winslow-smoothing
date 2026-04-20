function [p,t]=gri2pt(fn)

fid=fopen(fn,'r');

fgetl(fid);
a=fscanf(fid,'%d',3);
fgetl(fid);

nt=a(1);
np=a(2);

fgetl(fid);
t=zeros(nt,4);
for ii=1:nt
  t(ii,:)=fscanf(fid,'%d',4)';
end
fgetl(fid);

fgetl(fid);
p=zeros(np,3);
for ii=1:np
  p(ii,:)=fscanf(fid,'%f',3)';
end
fgetl(fid);

fgetl(fid);

fclose(fid);

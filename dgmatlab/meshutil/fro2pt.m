function [p,t]=fro2pt(fn)

fid=fopen(fn,'r');

a=fscanf(fid,'%d',6);

nt=a(1);
np=a(2);

p=zeros(np,3);
for ii=1:np
  foo=fscanf(fid,'%d',1);
  p(ii,:)=fscanf(fid,'%f',3)';
end

t=zeros(nt,4);
for ii=1:nt
  foo=fscanf(fid,'%d',1);
  t(ii,:)=fscanf(fid,'%d',4)';
end

fclose(fid);

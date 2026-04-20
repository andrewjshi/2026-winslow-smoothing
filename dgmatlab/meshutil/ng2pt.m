function [p,t,e,tnbr,enbr]=ng2pt(fname)
%NG2PT NETGEN mesh to P,T,E format
%   Currently assuming 3-D second-order elements
    
fid=fopen(fname,'r');
if fid==-1
    error('Can''t open file');
end

readuntil(fid,'surfaceelements');
ne=fscanf(fid,'%d',1);
e=zeros(ne,6);
enbr=zeros(ne,1);
for ii=1:ne
  enbr(ii)=fscanf(fid,'%d',1);
  foo=fscanf(fid,'%d',4);
  e(ii,:)=fscanf(fid,'%d',6)';
  fgetl(fid);
end

readuntil(fid,'volumeelements');
nt=fscanf(fid,'%d',1);
t=zeros(nt,10);
tnbr=zeros(nt,1);
for ii=1:nt
  tnbr(ii)=fscanf(fid,'%d',1);
  foo=fscanf(fid,'%d',1);
  t(ii,:)=fscanf(fid,'%d',10)';
  fgetl(fid);
end

readuntil(fid,'points');
np=fscanf(fid,'%d',1);
p=zeros(np,3);
for ii=1:np
  p(ii,:)=fscanf(fid,'%f',3)';
  fgetl(fid);
end

fclose(fid);

function readuntil(fid,str)
    
while ~feof(fid)
    fline=fgetl(fid);
    if ~isempty(fline) & ~isempty(strmatch(str, fline))
        break;
    end
end

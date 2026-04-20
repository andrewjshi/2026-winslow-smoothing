function fwritestruct(strct,fname,flds)

%if ~all(ismember(flds,fieldnames(strct)))
%  warning('Not all field names available in structure.');
%end

fid=fopen(fname,'w');
if fid==-1
  error('Can''t open file.');
end

for ii=1:numel(flds)
  if any(ismember(flds{ii},fieldnames(strct)))
    dat=getfield(strct,flds{ii});
  else
    dat=zeros(0,0);
  end
  fwritearray(fid,dat);
end

st=fclose(fid);
if st==-1
  error('Can''t close file.');
end

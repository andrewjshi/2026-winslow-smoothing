function strct=freadstruct(fname,flds)

fid=fopen(fname,'r');
if fid==-1
  error('Can''t open file.');
end

strct=struct;
for ii=1:numel(flds)
  strct=setfield(strct,flds{ii},freadarray(fid));
end

st=fclose(fid);
if st==-1
  error('Can''t close file.');
end

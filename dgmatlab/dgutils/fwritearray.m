function fwritearray(fid,a)

sz=size(a);
dim=numel(sz);
cl=class(a);

datatypes={'int32',0;
           'double',1;
           'logical',2};
ix=find(ismember(datatypes(:,1),cl));
if length(ix)~=1
  error('Datatype not supported.');
end
intcl=datatypes{ix,2};

isfilename=ischar(fid);
if isfilename
  fid=fopen(fid,'w');
end

fwrite(fid,[dim,sz,intcl],'int32');
fwrite(fid,a,cl);

if isfilename
  fclose(fid);
end

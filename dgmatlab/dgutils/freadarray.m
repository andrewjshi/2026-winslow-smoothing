function a=freadarray(fid)

isfilename=ischar(fid);
if isfilename
  fid=fopen(fid,'r');
end

dim=fread(fid,1,'int32');
sz=fread(fid,dim,'int32');
intcl=fread(fid,1,'int32');

datatypes={'int32',0;
           'double',1;
           'logical',2};
ix=find(ismember([datatypes{:,2}],intcl));
if length(ix)~=1
  error('Datatype not supported.');
end
cl=datatypes{ix,1};
mlcl=['*',cl];

szvec=sz(:)';
if numel(szvec)==1
    szvec(2)=1;
end
a=reshape(fread(fid,prod(sz),mlcl),szvec);

if isfilename
  fclose(fid);
end

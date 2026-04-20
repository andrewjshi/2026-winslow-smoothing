function createc(ccs,fname,mkheader,comments)

if nargin<3, mkheader=false; end
if nargin<4, comments={}; end

nccs=numel(ccs);
if numel(mkheader)==1
  mkheader=mkheader(ones(1,nccs));
end

hhs=cell(1,nccs);
funcs=cell(1,nccs);
for i=1:nccs
  cc=ccs{i};
  hh=extractheader(cc);
  cc=removemath(cc);
  func=extractfunc(cc);
  
  ccs{i}=cc;
  hhs{i}=hh;
  funcs{i}=func;
end

[foo,basefname,foo]=fileparts(fname);
headername=upper(basefname);
dotix=find(fname=='.');
hname=[fname(1:dotix(end)),'h'];

fid=fopen(fname,'w');
fprintf(fid,'#include <math.h>\n');
if any(mkheader)
  fprintf(fid,'#include "%s"\n',hname);
end
for i=1:nccs
  fprintf(fid,'\n');
  stars=repmat('*',1,78);
  fprintf(fid,'/%s/\n\n',stars);
  if nargin>=4 & ~isempty(comments{i})
    fprintf(fid,'/* %s: %s */\n',funcs{i},comments{i});
  end
  if ~any(mkheader) | mkheader(i), staticstr=''; else staticstr='static '; end
  fprintf(fid,'%s%s',staticstr,ccs{i});
end
fclose(fid);

if any(mkheader)
  fid=fopen(hname,'w');
  fprintf(fid,'#ifndef __%s\n',headername);
  fprintf(fid,'#define __%s\n',headername);
  for i=1:nccs
    if mkheader(i)
      fprintf(fid,'\n');
      if nargin>=4 & ~isempty(comments{i})
        fprintf(fid,'/* %s: %s */\n',funcs{i},comments{i});
      end
      fprintf(fid,'%s',hhs{i});
    end
  end
  fprintf(fid,'\n#endif\n');
  fclose(fid);
end

function cc=removemath(cc)

k1=findstr('void',cc);
if numel(k1)~=1 
  error('Assertion.');
end
cc=cc(k1:end);

function hh=extractheader(cc)

k1=findstr('{',cc);
k2=findstr('}',cc);
if numel(k1)==0 | numel(k2)==0
  error('Assertion.');
end

cc(k1(1):k2(end))=[];

k1=findstr('void',cc);
k2=findstr(')',cc);
if numel(k1)~=1 | numel(k2)~=1
  error('Assertion.');
end
hh=sprintf('%s;\n',cc(k1:k2));

function func=extractfunc(cc)

k1=findstr('void ',cc);
k2=findstr('(',cc);
if numel(k1)~=1 | numel(k2)<1
  error('Assertion.');
end
func=cc(k1+5:k2(1)-1);

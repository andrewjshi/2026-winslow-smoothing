function [ccode,ccode_deriv,maplein,mapleout]=autoc(in,pars,def,out,funcname,opts)
%AUTOC Automatic code generation, evaluation + derivatives
%
%  [CC,CC_DERIV]=AUTOC(IN,PARS,DEF,OUT,FUNCNAME,OPTSREMOVES,LANG,JACPERM)
%
%  OPTS.
%     REMOVES={REM_IN,REM_DEF,REM_OUT}
%     LANG = 'C'
%     JACPERM = {}
%     DOJAC = TRUE
%     DOMAPLE = TRUE
%
%  Example: 
%
%  [cc,cc_deriv]=autoc({'U',{'a','b'}},{},{'foo=a*sin(b)'}, ...
%                      {'F',{'foo'}},'calcfoo')

if nargin<6, opts=struct; end
if ~isfield(opts,'removes')
  opts.removes={cell(1,size(in,1)),cell(1,size(pars,1)),cell(1,size(out,1))};
end
if ~isfield(opts,'lang'), opts.lang='C'; end
if ~isfield(opts,'jacperm'), opts.jacperm={}; end
if ~isfield(opts,'dojac'), opts.dojac=true; end
if ~isfield(opts,'domaple'), opts.domaple=true; end

if nargout>=1
  [ccode,cmaplein,cmapleout]=subautoc(in,pars,def,out,funcname,opts,false);
  maplein{1}=cmaplein;
  mapleout{1}=cmapleout;
end
if nargout>=2
  if opts.dojac
    [ccode_deriv,cmaplein,cmapleout]=subautoc(in,pars,def,out,[funcname,'_deriv'],opts,true);
    maplein{2}=cmaplein;
    mapleout{2}=cmapleout;
  else
    ccode_deriv='';
  end
end

function [ccode,maplein,mapleout]=subautoc(in,pars,def,out,funcname,opts,jacobian)

predef={};
for i=1:numel(opts.removes{1})
  crem=opts.removes{1}{i};
  for j=1:numel(crem)
    predef{end+1}=[in{i,2}{crem(j)},'=0'];
  end
  in{i,2}(crem)=[];
end
for i=1:numel(opts.removes{2})
  crem=opts.removes{2}{i};
  for j=1:numel(crem)
    predef{end+1}=[pars{i,2}{crem(j)},'=0'];
  end
  pars{i,2}(crem)=[];
end
def={predef{:},def{:}}';

for i=1:numel(opts.removes{3})
  crem=opts.removes{3}{i};
  out{i,2}(crem)=[];
end

empty=sprintf('\n');

str=sprintf( ...
    ['with(linalg):\n', ...
     'with(codegen):\n', ...
     'with(CodeGeneration):\n', ...
     '\n']);

str=[str,sprintf( ...
    '`diff/abs` := proc(g,x) diff(g,x)*signum(g) end proc:\n')];

if isempty(pars), pars=cell(0,2); end

invar=in(:,1);
outvar=out(:,1);
parsvar=pars(:,1);
insize=varsz(in);
outsize=varsz(out);
parssize=varsz(pars);
if jacobian
  jacvar=mkjacvar(invar,outvar);
  jacsize=jacsz(insize,outsize);
end

str=[str,declare(in),empty];
str=[str,declare(pars),empty];

str=[str,def2maple(def),empty];

str=[str,mkcomputefunc(out),empty];

if jacobian
  str=[str,mkcomputejacs(invar,outvar,opts.jacperm),empty];
end

str=[str,mkprocs(invar,parsvar,outvar),empty];
if jacobian
  str=[str,mkjacprocs(jacvar),empty];
end

allvar=outvar;
if jacobian
  allvar={allvar{:},jacvar{:}};
end
allvarstr=sprintf('%s,',allvar{:});
str=[str,sprintf('%s:=makevoid(joinprocs([%s])):\n',funcname,allvarstr(1:end-1)),empty];

str=[str,mkdeclare(funcname,invar,insize)];
str=[str,mkdeclare(funcname,parsvar,parssize)];
str=[str,mkdeclare(funcname,outvar,outsize)];
if jacobian
  str=[str,mkdeclare(funcname,jacvar,jacsize)];
end
str=[str,empty];

tmpc=tempname;

str=[str,sprintf('%s(%s,optimize,deducetypes=false,output="%s"):\n', ...
                 opts.lang,funcname,tmpc)];

if opts.domaple
  tmpmaple=tempname;
  fid=fopen(tmpmaple,'w');
  fprintf(fid,'%s',str);
  fclose(fid);
  [status,mapleout]=system(sprintf('maple < %s',tmpmaple));
  if status~=0
    disp(mapleout);
    error('Error running maple, see output above.');
  end
  delete(tmpmaple);

  [foo,ccode]=system(sprintf('cat %s',tmpc));
  delete(tmpc);
else
  mapleout='';
  ccode='';
end
maplein=str;

function str=mkcomputefunc(out)

str='';
for i=1:size(out,1)
  v=sprintf('%s,',out{i,2}{:}); v(end)=[];
  str=[str,sprintf('%s:=evalm([%s]):\n',out{i,1},v)];
end

function str=mkdeclare(funcname,var,size)

str='';
for i=1:numel(var)
  str=[str,sprintf('%s:=declare(%s::''array''(1..%d,numeric),%s):\n', ...
                   funcname,var{i},size(i),funcname)];
end

function sz=varsz(var)

sz=zeros(size(var,1),1);
for i=1:size(var,1)
  sz(i)=length(var{i,2});
end
  
function sz=jacsz(insize,outsize)

sz=zeros(length(insize),length(outsize));
for i=1:length(insize)
  for j=1:length(outsize)
    sz(i,j)=insize(i)*outsize(j);
  end
end
    
function jacvar=mkjacvar(invar,outvar)

jacvar=cell(length(invar),length(outvar));
for i=1:length(invar)
  for j=1:length(outvar)
    jacvar{i,j}=sprintf('d%sd%s',outvar{j},invar{i});
  end
end

function str=declare(var)

str='';
for i=1:size(var,1)
  str=[str,sprintf('%s:=vector(%d):\n',var{i,1},length(var{i,2}))];
end
for i=1:size(var,1)
  for j=1:length(var{i,2});
    str=[str,sprintf('%s:=%s[%d]: ',var{i,2}{j},var{i,1},j)];
  end
  str=[str,sprintf('\n')];
end

function str=def2maple(def)

str='';
for i=1:size(def,1)
  cmd=def{i};
  str=[str,sprintf('%s:\n',strrep(cmd,'=',':='))];
end

function str=mkcomputejacs(invar,outvar,jacperm)

str='';
f=outvar;
u=invar;
for i=1:length(u)
  for j=1:length(f)
    str=[str,sprintf('d%sd%s:=jacobian(%s,%s):\n',f{j},u{i},f{j},u{i})];
  end
end
for i=1:length(u)
  for j=1:length(f)
    str=[str,sprintf('d%sd%s:=convert(evalm(transpose(d%sd%s)),vector):\n',f{j},u{i},f{j},u{i})];
  end
end
if ~isempty(jacperm)
  for i=1:length(u)
    for j=1:length(f)
      order=perm2order(jacperm{j,i}{:});
      orderstr=sprintf('%d,',order); orderstr(end)=[];
      str=[str,sprintf('d%sd%s:=vector([seq(d%sd%s[i],i=[%s])]):\n',f{j},u{i},f{j},u{i},orderstr)];
    end
  end
end

function order=perm2order(sz,prm)

order=reshape(permute(reshape(1:prod(sz),sz),prm),[],1);

function str=mkprocs(invar,parsvar,outvar)

str='';
invars=sprintf('%s,',invar{:});
parsvars=sprintf('%s,',parsvar{:});
for i=1:length(outvar)
  v=outvar{i};
  if i==1
    dep=sprintf('[%s%s%s]',invars,parsvars,v);
  else
    dep=v;
  end
  str=[str,sprintf('%s:=makevoid(makeproc(%s,%s)):\n',v,v,dep)];
end

function str=mkjacprocs(jacvar)

str='';
for i=1:numel(jacvar)
  v=jacvar{i};
  str=[str,sprintf('%s:=makevoid(makeproc(%s,%s)):\n',v,v,v)];
end

function dgexport(format,fname,msh,u,nref,navstok,solutiontime)
%DGEXPORT Export DG data to Tecplot.
%   DGEXPORT(FORMAT,FNAME,MSH,U,NREF,NAVSTOK,SOLTIME)
%
% DG BINARY FILE FORMAT
% ---------------------
%
%  D                         1 int32    number of space dimensions
%  np                        1 int32    number of node points
%  nt                        1 int32    number of volume simplex elements
%  nb                        1 int32    number of boundary simplex elements
%  nc                        1 int32    number of visualization components
%  x1,y1,z1,...              D single   node point 1 coordinates
%  c11,c21,c31,...          nc single   node point 1 component values
%  ...                                  repeat for np node points
%  n11,n21,n31,n41,...     D+1 int32    volume element 1 node indices
%  ...                                  repeat for nt volume elements
%  b11,b21,b31,bnbr1,...   D+1 int32    boundary element 1 indices + nbr
%  ...                                  repeat for nb boundary elements
%
%   Example: Export Naiver-Stokes problem to V3, no extra refinements
%      dgexport('v3','model1.v3',msh,u,0,true);

if nargin<5, nref=0; end
if nargin<6, navstok=false; end
if nargin<7, solutiontime=[]; end

np=size(msh.p,2);
nt=size(msh.t,2);
ns=size(msh.s,1);
dim=size(msh.p,1);

porder=double(msh.porder);

if porder==0
  [ss,tt]=lagrangepnts(1,dim);
  pp=zeros(dim+1,dim,nt);
  for i=1:dim
    pp(:,i,:)=reshape(msh.p(i,msh.t+1),dim+1,1,nt);
  end
  u=repmat(u,[dim+1,1,1]);
else
  ss=msh.s;
  tt=msh.tlocal;
  pp=msh.p1;
end

if min(tt(:)) == 0 % Support zero-based format
    tt = tt + 1;
end

if size(u,1)==1
  u=repmat(u,[ns,1,1]);
end

if nref>0
  if dim==2
    A0=pmonomial(ss(:,1:2),porder);
    [ss,tt]=uniref(ss,tt,nref);
    A=pmonomial(ss(:,1:2),porder)/A0;
    
    nss=size(ss,1);
    sz=size(pp); if length(sz)==2, sz=[sz,1]; end
    pp=reshape(A*reshape(pp,ns,sz(2)*sz(3)),[nss,sz(2),sz(3)]);
    sz=size(u); if length(sz)==2, sz=[sz,1]; end
    u=reshape(A*reshape(u,ns,sz(2)*sz(3)),[nss,sz(2),sz(3)]);
  elseif dim==3
    A0=pmonomial(ss(:,1:3),porder);
    [ss,tt]=lagrangepnts(porder*2^nref,3);
    A=pmonomial(ss(:,1:3),porder)/A0;
    
    nss=size(ss,1);
    sz=size(pp); if length(sz)==2, sz=[sz,1]; end
    pp=reshape(A*reshape(pp,ns,sz(2)*sz(3)),[nss,sz(2),sz(3)]);
    sz=size(u); if length(sz)==2, sz=[sz,1]; end
    u=reshape(A*reshape(u,ns,sz(2)*sz(3)),[nss,sz(2),sz(3)]);
  end
end

nss=size(ss,1);
p1vis=reshape(permute(pp,[1,3,2]),[nss*nt,dim]);
uvis=reshape(permute(u,[1,3,2]),[nss*nt,size(u,2)]);
tvis=kron(ones(nt,1,'int32'),tt)+kron(nss*int32(0:nt-1)',0*tt+1);

[p1vis,tvis,p1ix]=fixmesh(p1vis,tvis,1e-5);
uvis=uvis(p1ix,:);

spacevars={'x','y','z'};
fieldvars={'u1','u2','u3','u4','u5','u6','u7','u8','u9'};
extravars={};
extravis=zeros(nss*nt,0);
if navstok
  if dim==2
    fieldvars(1:4)={'r','ru','rv','rE'};
  else
    fieldvars(1:5)={'r','ru','rv','rw','rE'};
  end
  
  if navstok==2
    extravars={'u','v','w','p','M','s','uabs'};
  
    extra_u=dgeval(u,'u',dim);
    extra_v=dgeval(u,'v',dim);
    extra_w=dgeval(u,'w',dim);
    extra_p=dgeval(u,'p',dim);
    extra_M=dgeval(u,'M',dim);
    extra_s=dgeval(u,'s',dim);
    extra_uabs=dgeval(u,'uabs',dim);
    extra=cat(2,extra_u,extra_v,extra_w,extra_p,extra_M,extra_s,extra_uabs);
    extravis=reshape(permute(extra,[1,3,2]),[nss*nt,size(extra,2)]);
  end
end
extravis=extravis(p1ix,:);

if max(abs(imag(uvis)))~=0, warning('Imaginary data'); end
uvis=real(uvis);
if max(abs(imag(extravis)))~=0, warning('Imaginary data'); end
extravis=real(extravis);

nc=size(uvis,2)+size(extravis,2);
dat=[p1vis,uvis,extravis];

np=size(p1vis,1);
nt=size(tvis,1);

if dim==2
  nb=0;
  tri=[];
  bndnbrs=[];
elseif dim==3
  tri=surftri(p1vis,tvis);
  nb=size(tri,1);

  if isfield(msh,'bndexpr') & ~isempty(msh.bndexpr)
    bndnbrs=zeros(nb,1);
    for ib=1:nb
      p=p1vis(tri(ib,:),:);
      found=false;
      for jj=1:length(msh.bndexpr)
        if eval(msh.bndexpr{jj})
          found=true;
          bnd=jj;
          break;
        end
      end
      if ~found, error('Strange boundary.'); end
      bndnbrs(ib)=bnd;  
    end
  else
    bndnbrs=ones(nb,1);
  end
end

switch format
 case 'v3'
  fid=fopen(fname,'w');
  fwrite(fid,int32([dim,np,nt,nb,nc]),'int32');
  %fwrite(fid,single(dat)','single');
  fwrite(fid,single(p1vis)','single');
  fwrite(fid,single(uvis)','single');
  fwrite(fid,tvis','int32');
  fwrite(fid,int32([tri,bndnbrs])','int32');
  fclose(fid);
 case 'tecplot'
  fid=fopen(fname,'w');
  fprintf(fid,'title="DG"\n');
  fprintf(fid,'variables=');
  allvars={spacevars{1:dim},fieldvars{1:size(u,2)},extravars{:}};
  for i=1:numel(allvars)
    if i>1, fprintf(fid,','); end
    fprintf(fid,'"%s"',allvars{i});
  end
  fprintf(fid,'\n');
  elementnames={'','fetriangle','fetetrahedron'};
  fprintf(fid,'zone T="Volume" n=%d, e=%d, datapacking=point, zonetype=%s\n', ...
          size(p1vis,1),size(tvis,1),elementnames{dim});
  if ~isempty(solutiontime)
      fprintf(fid,'SolutionTime=%.6g\n',solutiontime);
  end
  dat=[p1vis,uvis,extravis];
  fmtstr = repmat('%.6g,', 1, size(dat,2));
  fmtstr = [fmtstr(1:end-1),'\n'];
  fprintf(fid, fmtstr, dat');
  fprintf(fid, '%d,%d,%d,%d\n', tvis');
  if dim==3
    Dstr='';
    for i=1:numel(allvars)
      Dstr=[Dstr,int2str(i),','];
    end
    Dstr(end)=[];
    
    for ibnd=1:max(bndnbrs)
      zonename=sprintf('Surface %d',ibnd);
      ctri=tri(bndnbrs==ibnd,:);
      fprintf(fid,['zone T="%s" n=%d, e=%d, datapacking=point, ', ...
                   'zonetype=fetriangle, D=(%s)\n'], ...
              zonename,size(p1vis,1),size(ctri,1),Dstr);
      fprintf(fid, '%d,%d,%d\n', ctri');
    end
  end
  fclose(fid);
 case 'vtk'
  fid=fopen(fname,'w');
  fprintf(fid,'# vtk DataFile Version 2.0\n');
  fprintf(fid,'DG\n');
  fprintf(fid,'ASCII\n');
  fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
  fprintf(fid,'POINTS %d float\n',size(p1vis,1));
  fclose(fid);
  xyz = p1vis;
  if dim == 2
      xyz = [xyz, 0*xyz(:,1)];
  end
  dlmwrite(fname, xyz, '-append', 'delimiter',' ');

  fid=fopen(fname,'a');
  fprintf(fid,'\nCELLS %d %d\n',size(tvis,1),(1 + size(tvis,2)) * size(tvis,1));
  fclose(fid);
  dlmwrite(fname,[0*tvis(:,1)+size(tvis,2),tvis-1],'-append','precision','%d','delimiter',' ');

  fid=fopen(fname,'a');
  fprintf(fid,'\nCELL_TYPES %d\n',size(tvis,1));
  fclose(fid);
  if dim == 2
      celltype = 5; % triangle
  elseif dim == 3
      celltype = 10; % tet
  end
  dlmwrite(fname, 0*tvis(:,1) + celltype, '-append', 'precision','%d', 'delimiter',' ');

  fid=fopen(fname,'a');
  fprintf(fid,'\nPOINT_DATA %d\n',size(uvis,1));
  fclose(fid);
    
  if ~navstok % General case, just dump all the solution components
      mk_vtk_data(fname, 'SCALARS', 'scalar_data', uvis)
  else % Navier-Stokes
      mk_vtk_data(fname, 'SCALARS', 'Density_r', uvis(:,1));
      if size(uvis,2) >= dim+2
          mk_vtk_data(fname, 'SCALARS', 'Energy_rE', uvis(:,dim+2));
      end
      mk_vtk_data(fname, 'VECTORS', 'Momentum_ru', uvis(:,2:dim+1));
  end      
  otherwise
    error('Unknown export format.');
end

end

function mk_vtk_data(fname, typename, uname, u)
    fid=fopen(fname,'a');
    if all(typename=='VECTORS')
        fprintf(fid,'%s %s float\n', typename, uname);
        if size(u,2) == 2
            u = [u, 0*u(:,1)];
        end
    else
        fprintf(fid,'%s %s float %d\n', typename, uname, size(u,2));
        fprintf(fid,'LOOKUP_TABLE default\n');
    end
    fclose(fid);
    dlmwrite(fname,u,'-append','delimiter',' ');
end

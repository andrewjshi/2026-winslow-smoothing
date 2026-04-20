function msh=mshpartition(msh,nregion)
%MSHPARTITION  Partition DG mesh.
%
%    Syntax: msh=mshpartition(msh,nregion)

t2t=msh.t2t'+1;
A=dualadjacency(t2t);
nt=size(t2t,1);

if isscalar(nregion)
  region=metismex('PartGraphRecursive',A,nregion)+1;
  if all(region==2) % Seems like a bug in metis
    region(:)=1;
  end
else % assume partition input
  region=nregion(:)';
  nregion=max(region(:));
end
local=zeros(1,nt);

partsz=zeros(3,nregion);
BINs=cell(1,nregion);
iBIs=cell(1,nregion);
cellBINs=cell(1,nregion);
lt2ts=cell(1,nregion);
comm=zeros(5,0,'int32');

for i=1:nregion
  ix=find(region==i);
  ct2t=t2t(ix,:);
  posix=find(ct2t>=1);
  neighbor_region=0*ct2t;
  neighbor_region(posix)=region(ct2t(posix));
  ixB=find(any(neighbor_region>=1 & neighbor_region~=i,2))';
  ixI=setdiff(1:length(ix),ixB);
  B=ix(ixB);
  I=ix(ixI);
  N=setdiff(ct2t(:)',[B,I,min(t2t(:)):0]);
  
  Bt2t=ct2t(ixB,:);
  Nix=ismember(Bt2t,N);
  Bexp=repmat(B(:),1,size(t2t,2));
  
  comm_from=Bexp(Nix);
  comm_to=Bt2t(Nix);
  comm_to_region=region(comm_to)';
  unique_comm_to_region=unique(comm_to_region);
  X=sparse(comm_from,comm_to_region,1);
  X=~~full(X(B,unique_comm_to_region));
  pp=Boptorder(X);
  B=B(pp);
  X=X(pp,:);
  
  cproc=find(sum(X,1)>0);
  for j=cproc
    cncomm=reshape(find(diff([0;X(:,j);0])),2,[]);
    cncomm(2,:)=cncomm(2,:)-cncomm(1,:);
    for k=1:size(cncomm,2)
      comm(:,end+1)=int32([i-1; unique_comm_to_region(j)-1;
                          cncomm(1,k)-1; 0; cncomm(2,k)]);
    end
  end
  
  cellBINs{i}={B,I,N};
end

for i=1:nregion
  B=cellBINs{i}{1};
  I=cellBINs{i}{2};
  N=cellBINs{i}{3};
  Nregions=region(N);
  [foo,sortNix]=sort(Nregions);
  N=N(sortNix);


  newN=[];
  if ~isempty(Nregions)
    for j=unique(Nregions)
      Bneighbor=cellBINs{j}{1};
      myBneighbors=Bneighbor(ismember(Bneighbor,N));
      newN=[newN,myBneighbors];
      allk=find(comm(1,:)==j-1 & comm(2,:)==i-1);
      for k=allk
        comm(4,k)=int32(find(newN==Bneighbor(comm(3,k)+1))-1);
      end
    end
    N=newN;  % newN seems to be equal to N anyway because of sorting
  end

  local([B,I])=1:(length(B)+length(I));
  glob=[B,I,N];
  
  map=-999999*ones(1,nt);
  map(glob)=1:length(glob);
%  ct2t=ct2t([ixB,ixI],:);
  ct2t=t2t(glob,:);
  posix=find(ct2t>=1);
  lt2t=ct2t;
  lt2t(posix)=map(ct2t(posix));

  partsz(:,i)=cumsum([length(B);length(I);length(N)]);
  BINs{i}=int32(glob)-1;
  lt2ts{i}=int32(lt2t)'-1;
  [foo,iBI]=sort([B,I]);
  iBIs{i}=int32(iBI)-1;
end

msh.tpart=int32([region;local])-1;
msh.partsz=int32([partsz;[0;0],cumsum(partsz(2:3,1:nregion-1),2)]);
msh.BINs=[BINs{:}];
msh.iBIs=[iBIs{:}];
msh.lt2ts=[lt2ts{:}];
msh.comm=comm;

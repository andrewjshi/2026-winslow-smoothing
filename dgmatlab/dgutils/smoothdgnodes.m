function msh=smoothdgnodes(msh)

if ~(isfield(msh,'fd') & ~isempty(msh.fd))
  error('No distance function available.');
end

porder=msh.porder;
dim=size(msh.p,1);

fd=msh.fd;
if isfield(msh,'fdargs')
  fdargs=msh.fdargs;
else
  fdargs={};
end

if 0 & all(msh.tcurved)
  % Special technique: Project boundary nodes only to level set
  % Elliptic smoothing for all interior nodes

  newp1=msh.p1;
  tix=find(any(msh.t2t<0,1));
  for it=tix
    p=msh.p1(:,:,it);
    isbnd=msh.ncurved(msh.t(:,it)+1);
    isbnd=repmat(isbnd,size(msh.s,1),1);
    bix=all(isbnd | msh.s<1e-6,2);
    
    if any(bix)
      cix=sum(msh.s>1-1e-6,2)==1;
      cix=cix|~bix;
      
      ns=size(msh.s,1);
      S=speye(ns,ns);
      if dim==2 % Some smoothing along boundary edge
        ee=boundedges(msh.s(:,1:2),double(msh.tlocal+1));
        ee=ee(all(bix(ee),2),:);
        degree=full(sparse(ee,1,1,ns,1));
        S=sparse(ee,fliplr(ee),1./degree(ee),ns,ns);
        S(cix,:)=0;
        S(cix,cix)=speye(sum(cix),sum(cix));
      end
      
      for iii=1:5
        %      p(bix,:)=bndproj(p(bix,:),msh.tlocal+1,fd,fdargs{:});
        p=S*p;
        p(bix,:)=pproj(p(bix,:),0,fd,fdargs{:});
      end

      newp1(:,:,it)=p;
    end
  end

  msh=elliptic_smoothing(msh,newp1,1);
  
  % Special technique: Project to level sets
  % To do: Handle only part of the elements this way
%  for it=1:size(msh.t,2)
%    p1=msh.p1(:,:,it);
%    p=msh.p(:,msh.t(:,it)+1)';
%    d=feval(fd,p,fdargs{:});
%    d0=msh.s*d;
%    for iii=1:5
%      p1=pproj(p1,d0,fd,fdargs{:});
%    end
%    msh.p1(:,:,it)=p1;
%  end
else
  % Typical boundaries: Project boundary nodes only to level set
  % Smooth remaining nodes
  tix=find(msh.tcurved);
  for it=tix
    p=msh.p1(:,:,it);
    
    if msh.eltype==t_simplex
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Simplex elements
        
        isbnd=msh.ncurved(msh.t(:,it)+1);
        isbnd=repmat(isbnd,size(msh.s,1),1);
        bix=all(isbnd | msh.s<1e-6,2);
        cix=sum(msh.s>1-1e-6,2)==1;
        cix=cix|~bix;
        
        ns=size(msh.s,1);
        S=speye(ns,ns);
        if dim==2 % Some smoothing along boundary edge
            ee=boundedges(msh.s(:,1:2),double(msh.tlocal+1));
            ee=ee(all(bix(ee),2),:);
            degree=full(sparse(ee,1,1,ns,1));
            S=sparse(ee,fliplr(ee),1./degree(ee),ns,ns);
            S(cix,:)=0;
            S(cix,cix)=speye(sum(cix),sum(cix));
        end
        
        if any(bix)
            for iii=1:3
                %      p(bix,:)=bndproj(p(bix,:),msh.tlocal+1,fd,fdargs{:});
                p=S*p;
                p(bix,:)=pproj(p(bix,:),0,fd,fdargs{:});
            end
        end
        
        % Smooth
        % Smoothing currently only implemented for 2-D
        % To Do: In 3-D, first smooth faces, then interior
        
        t=double(msh.tlocal+1);
        if dim==2
            bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];
        else
            bars=[t(:,[1,2]);t(:,[1,3]);t(:,[1,4]);t(:,[2,3]);t(:,[2,4]);t(:,[3,4])];
        end
        bars=unique(sort(bars,2),'rows');
        
        if dim==2
            nsm=10*porder;
            iix=find(all(msh.s>1e-6,2));
            for ism=1:nsm
                ps=p;
                for idim=1:dim
                    cp=p(:,idim);
                    ps(:,idim)=full(sparse(bars,1,fliplr(cp(bars)),size(p,1),1)./ ...
                                    sparse(bars,1,1,size(p,1),1));
                end
                p(iix,:)=ps(iix,:);
            end
        end
    elseif msh.eltype==t_block
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Block elements
        
        if dim==2
            p=reshape(p,[porder+1,porder+1,dim]);
            p0=p;
            
            e1=permute(p(:,1,:),[1,3,2]);
            e2=permute(p(end,:,:),[2,3,1]);
            e3=permute(p(:,end,:),[1,3,2]);
            e4=permute(p(1,:,:),[2,3,1]);
            es={e1,e2,e3,e4};

            xy=msh.s*[0,0;1,0;0,1;1,1];
            s=xy(xy(:,2)<1e-6,1);
            kmap=[3,2,4,1];
            for k=1:4
                if msh.ecurved(kmap(k),it)
                    for i=1:5
                        ss=sqrt(sum(diff(es{k},[],1).^2,2));
                        es1=es{k};
                        es{k}=interp1([0;cumsum(ss)],es{k},s*sum(ss),'linear','extrap');
                        es2=es{k};
                        es{k}=pproj(es{k},0,fd,fdargs{:});
                        es3=es{k};
                    end
                end
            end
            
            % TFI
            x1=es{1}(:,1); y1=es{1}(:,2);
            x2=es{3}(:,1); y2=es{3}(:,2);
            x3=es{4}(:,1); y3=es{4}(:,2);
            x4=es{2}(:,1); y4=es{2}(:,2);
            
            X1=(1-s)*x1'+s*x2';
            Y1=(1-s)*y1'+s*y2';
            X2=x3*(1-s')+x4*s';
            Y2=y3*(1-s')+y4*s';
            X12=(1-s)*(1-s')*X1(1,1)+s*(1-s')*X1(end,1)+(1-s)*s'*X1(1,end)+s*s'*X1(end,end);
            Y12=(1-s)*(1-s')*Y1(1,1)+s*(1-s')*Y1(end,1)+(1-s)*s'*Y1(1,end)+s*s'*Y1(end,end);
            X=X1+X2-X12;
            Y=Y1+Y2-Y12;

            p(:,:,1)=rot90(fliplr(X));
            p(:,:,2)=rot90(fliplr(Y));

            p=reshape(p,[(porder+1)^2,dim]);
        else
            error('not implemented');
        end
    end
        
    msh.p1(:,:,it)=p;
  end
end

function varargout=froedgeinterp(fro,nn,s,onlycurves)

varargout=cell(1,nargout);

s=s(:,1);
s=[s,1-s];

npc=numel(fro.pc);
for ipc=1:npc
  cpc=fro.pc{ipc};
  if any(cpc(:,1)==nn(1)) && any(cpc(:,1)==nn(2))
    ix1=find(cpc(:,1)==nn(1));
    ix2=find(cpc(:,1)==nn(2));
    u=s*cpc([ix1,ix2],2);
    [varargout{:}]=froceval(fro,ipc,u);
    return
  end
end

if nargin>=4 && onlycurves
  [varargout{:}]=deal(nan);
  return;
end

nps=numel(fro.ps);
for ips=1:nps
  cps=fro.ps{ips};
  if any(cps(:,1)==nn(1)) && any(cps(:,1)==nn(2))
    ix1=find(cps(:,1)==nn(1));
    ix2=find(cps(:,1)==nn(2));
    v=s*cps([ix1,ix2],2);
    w=s*cps([ix1,ix2],3);
    [varargout{:}]=froseval(fro,ips,[v,w]);
    return
  end
end

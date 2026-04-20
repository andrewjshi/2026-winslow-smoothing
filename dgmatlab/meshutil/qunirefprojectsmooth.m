function [p,q]=qunirefprojectsmooth(p,q,nref,nsmooth,fd,varargin)

for iref=1:nref
    [p,q]=quniref(p,q,1);
    bndix=unique(boundedges(p,[q(:,[1,2,3]);q(:,[2,4,3])]));
    p(bndix,:)=pproj(p(bndix,:),0,fd,varargin{:});
    p=quadsmooth(p,q,bndix,nsmooth);
end

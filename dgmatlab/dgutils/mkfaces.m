function faces=mkfaces(t,eltype)
    
if nargin<2, eltype=t_simplex; end
dim=nv2dim(size(t,2),eltype);

facemap=mkfacemap(dim,eltype);
faces=[];
for i=1:size(facemap,2)
    faces=[faces; t(:,facemap(:,i))];
end

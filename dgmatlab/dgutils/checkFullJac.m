if ~exist('it'), it=1; end
if ~exist('resname'), resname=@dgnavierstokes; end

[R0,DJ0,OJ0]=feval(resname,u,msh,data,phys);
J=mkjacobian(msh.t2t,data.egix,DJ0,OJ0);

n=size(R0,1)*size(R0,2);
nu=size(u,1)*size(u,2);

deltau=1e-6;
R0=feval(resname,u,msh,data,phys);
D0=zeros(prod(size(R0)),n);
D1=full(J(:,n*(it-1)+1:n*it));
for i1=1:size(R0,1)
    for i2=1:size(R0,2)
        u1=u;
        u1(i1,i2,it)=u1(i1,i2,it)+deltau;
        R1=feval(resname,u1,msh,data,phys);
        u2=u;
        u2(i1,i2,it)=u2(i1,i2,it)-deltau;
        R2=feval(resname,u2,msh,data,phys);
        cD=(R1-R2)/2/deltau;
        D0(:,i1+(i2-1)*size(R0,1))=cD(:);
    end
end
fprintf('%d %g\n',it,max(max(abs(D1-D0))));

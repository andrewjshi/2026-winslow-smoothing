function XH=block2hexes(X)
%BLOCK2HEXES  Subdivides a block into hexes.
%
%   XH=block2hexes(X)

dim = size(X,1);

switch dim
    case 1
        XH = [X(1:end-1);X(2:end)];
        XH = reshape(XH,1,2,[]);
    case 2
        nx = size(X,2); ny = size(X,3);
        XH = zeros(2,4,(nx-1)*(ny-1));
        ik = 1;
        for j=1:ny-1
            for i=1:nx-1
                XH(1:2,1:4,ik) = [
                    X(:,i,j  ), X(:,i+1,j  ), ...
                    X(:,i,j+1  ), X(:,i+1,j+1  ) ];
                ik = ik+1;
            end
        end
        XH=XH(:,[1,3,2,4],:); % weird orientation, to be consistent with tets
    case 3
        nx = size(X,2); ny = size(X,3); nz = size(X,4);
        XH = zeros(3,8,(nx-1)*(ny-1)*(nz-1));
        ik = 1;
        for k=1:nz-1
            for j=1:ny-1
                for i=1:nx-1
                    XH(1:3,1:8,ik) = [
                        X(:,i,j,k  ), X(:,i+1,j,k  ), ...
                        X(:,i,j+1,k  ), X(:,i+1,j+1,k  ), ...
                        X(:,i,j,k+1), X(:,i+1,j,k+1), ...
                        X(:,i,j+1,k+1), X(:,i+1,j+1,k+1) ];
                    ik = ik+1;
                end
            end
        end
end

function up = dginterp2d(msh,X0,u,verbose)

if nargin<4, verbose=true; end

eps = 0.1;               % Lenght margin (increase for robustness - should nto be larger than 0.2)
nref = 0;                % Number of refinements for initial guess  (increase for robustness)

dg = msh.p1;

np = size(X0,1);
ne = size(msh.t2n,2);

porder = double(msh.porder);
s = msh.s;
A0=pmonomial(s(:,1:2),porder);
ns=size(s,1);

pf = porder*2^nref;
[uc,vc]=ndgrid((0:pf)/pf,(0:pf)/pf);
ss=[uc(:),vc(:)];
ss=[ss,1-sum(ss,2)];
ss=ss(ss(:,3)>=0,:);

A=pmonomial(ss(:,1:2),porder)/A0;

up = zeros(np,size(u,2));

for ip = 1:np
    p = X0(ip,:);
    
    if 0
        % Original
        elm = zeros(ne,1);
        for ie = 1:ne         % Eliminate elements that are irrelevant
            epsx = eps*(max(dg(:,1,ie))-min(dg(:,1,ie)));
            epsy = eps*(max(dg(:,2,ie))-min(dg(:,2,ie)));
            elm(ie) = all(dg(:,1,ie)>p(:,1)+epsx) || all(dg(:,2,ie)>p(:,2)+epsy) ||  ...
                      all(dg(:,1,ie)<p(:,1)-epsx) || all(dg(:,2,ie)<p(:,2)-epsy) ;
        end
    elseif 0
        % Vectorized
        epsx=eps*(max(dg(:,1,:),[],1)-min(dg(:,1,:),[],1));
        epsy=eps*(max(dg(:,2,:),[],1)-min(dg(:,2,:),[],1));
        elm=all(dg(:,1,:)>p(:,1)+repmat(epsx,[ns,1,1]),1) | ...
            all(dg(:,2,:)>p(:,2)+repmat(epsy,[ns,1,1]),1) | ...
            all(dg(:,1,:)<p(:,1)-repmat(epsx,[ns,1,1]),1) | ...
            all(dg(:,2,:)<p(:,2)-repmat(epsy,[ns,1,1]),1);
    else
        % Mexed
        elm=dginterp_elim(dg,p);
    end
    list = find(elm == 0);
    dist = zeros(length(list),1);
   
    for i = 1:length(list)   % Refine element nref times and get distance to closest dgnode
        ie = list(i);
        xc = A*dg(:,:,ie);
        dis = (xc(:,1)-p(:,1)).^2 + (xc(:,2)-p(:,2)).^2;
        dist(i) = min(dis);
    end
    
    [dum,il] = sort(dist);   % Sort elements into list of increasing likelyhood
    list = list(il);

    found = false;
    for i = 1:length(list)   % Run trough list
        ie = list(i);
        xc = A*dg(:,:,ie);
        dis = (xc(:,1)-p(:,1)).^2 + (xc(:,2)-p(:,2)).^2;
        [d,ind] = min(dis);
        sc = ss(ind(1),1:2); % Closest dg node
        d = 1;
        iter = 0;
        while d > 1.e-10     % Newton iteration
           [B,Bx,By] = pmonomial(sc,porder);
           B = B/A0; Bx = Bx/A0; By = By/A0;;
           xc = B*dg(:,:,ie);
           J = [Bx*dg(:,:,ie); By*dg(:,:,ie)]';
           R = (p-xc)';
           dsc = J\R;
           sc = sc + dsc';
           d = norm(R);
           iter = iter + 1;
           if iter>10, break; end
        end
        up(ip,:) = B*u(:,:,ie);  % Interpolate unknown
        
        sc=[sc,1-sum(sc,2)];
        if iter < 11 & max(sc) < 1.00000001 & min(sc) > -0.0000001
            found = true;
            break;
        end   %Make sure we are inside element
    end

    if ~found
        up(ip,:) = nan;
        continue;
    end
    
    if verbose
        fprintf('Point %d - iter %d  - itry %d - scmax %f - scmin %f \n',ip,iter,i,max(sc),min(sc));
    end
end


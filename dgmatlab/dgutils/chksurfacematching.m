function chksurfacematching(msh1, bndnbr1, msh2, bndnbr2)
% CHKSURFACEMATCHING Validate the matching of the surfaces.
%   CHKSURFACEMATCHING(MSH1, BNDNBR1, MSH2, BNDNBR2) throws an
%   error if the faces of boundary number BNDNBR1 of mesh MSH1 do
%   not match with the faces of boundary number BNDNBR2 of mesh
%   MSH2.

dim = size(msh1.p,1);
eltype = msh1.eltype;
if eltype ~= t_simplex
  error('Not implemented.');
end
faces = mkfacemap(dim,eltype);
nf = size(msh1.t2t,1);
nv1 = size(faces,1);

p1 = msh1.p';
t1 = msh1.t'+1;
t2t1 = msh1.t2t';

p2 = msh2.p';
t2 = msh2.t'+1;

% Extract boundary triangulation.
et1 = zeros(0, nv1);
[j1, it1]=find(ismember(-msh1.t2t, bndnbr1));
for i=1:numel(it1)
  et1 = [et1; t1(it1(i), faces(:,j1(i)))];
end

et2 = zeros(0, nv1);
[j2, it2]=find(ismember(-msh2.t2t, bndnbr2));
for i=1:numel(it2)
  et2 = [et2; t2(it2(i), faces(:,j2(i)))];
end

if ~all(size(et1) == size(et2))
  error('Different number of boundary faces');
end

% Renumber and throw away unused points.
[ix,~,jx] = unique(et1(:));
p1 = p1(ix,:);
et1 = reshape(jx,size(et1));

[ix,~,jx] = unique(et2(:));
p2 = p2(ix,:);
et2 = reshape(jx,size(et2));

% Brute force map between p1 and p2:
d = sqrt(sum((repmat(reshape(p1, [size(p1,1),1,dim]),[1,size(p1,1),1]) ...
    -repmat(reshape(p2, [1,size(p2,1),dim]),[size(p2,1),1,1])).^2, 3));
[i1,i2] = find(d < 1e-8);
[~,j1] = sort(i1);
[~,j2] = sort(i2);

if (numel(i1) ~= size(p1,1))
  error('Points do not all align');
end
% At this point p1(i1,:) == p2(i2,:), at least up to tolerance.

% Find midpoints
pmid1 = p1(et1(:,1),:);
pmid2 = p2(et2(:,1),:);
for i=2:dim
    pmid1 = pmid1 + p1(et1(:,i),:);
    pmid2 = pmid2 + p2(et2(:,i),:);
end
pmid1 = pmid1/dim;
pmid2 = pmid2/dim;

d = sqrt(sum((repmat(reshape(pmid1, [size(pmid1,1),1,dim]),[1,size(pmid1,1),1]) ...
    -repmat(reshape(pmid2, [1,size(pmid2,1),dim]),[size(pmid2,1),1,1])).^2, 3));
[imid1,imid2] = find(d < 1e-8);

if (numel(imid1) ~= size(et2,1))
  error('Midpoints do not all align');
end

a1 = j1(et1(imid1,:));
a2 = j2(et2(imid2,:));

% STUPID MATLAB NEEDS ACTUAL VECTORS
% workaround breakage if the boundary contains only one face.
if size(a1,2) == 1, a1 = a1'; end
if size(a2,2) == 1, a2 = a2'; end

match = isequal(sort(a1,2), sort(a2,2));
if ~match
    error('Faces do not match');
end

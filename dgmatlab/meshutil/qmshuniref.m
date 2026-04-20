function msh1 = qmshuniref(msh, nref)
% Uniform high-order quad-mesh refinement

msh1 = msh;
for iref = 1:nref
    msh1 = qmshrefine(msh1, true(size(msh1.t')));
end

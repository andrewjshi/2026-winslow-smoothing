function out = ac_isentrop_inviscid(U, Fi, GAMMA, k_pr)
% INVISCID Fluxes
% out = ac_isentrop_inviscid(U, Fi, GAMMA, k_pr)

if (isnumeric(GAMMA))
  GAMMA = num2str(GAMMA);
end
if (isnumeric(k_pr))
  k_pr = num2str(k_pr);
end

out = { ...
'p=[k_pr]*[r]^([GAMMA])';
'[F]11=[ru]';
'[F]21=[ru]^2/[r]+p';
'[F]31=[ru]*[rv]/[r]';
'[F]41=[ru]*[rw]/[r]';
'[F]12=[rv]';
'[F]22=[ru]*[rv]/[r]';
'[F]32=[rv]^2/[r]+p';
'[F]42=[rv]*[rw]/[r]';
'[F]13=[rw]';
'[F]23=[ru]*[rw]/[r]';
'[F]33=[rv]*[rw]/[r]';
'[F]43=[rw]^2/[r]+p';
};

replacements = {
% U
'[r]',  [U '1'];
'[ru]', [U '2'];
'[rv]', [U '3'];
'[rw]', [U '4'];
% Fi
'[F]',  Fi;
% GAMMA
'[GAMMA]', GAMMA;
% k_pr
'[k_pr]', k_pr;
};

for r=replacements'
  out = strrep(out,r{1},r{2});
end


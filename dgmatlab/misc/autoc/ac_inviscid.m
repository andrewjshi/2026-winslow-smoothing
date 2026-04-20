function out = ac_inviscid(U, Fi, GAMMA)
% INVISCID Fluxes
% out = ac_inviscid(U, Fi, GAMMA)

if (isnumeric(GAMMA))
  GAMMA = num2str(GAMMA);
end

out = { ...
'p=([GAMMA]-1)*([rE]-1/2*([ru]^2+[rv]^2+[rw]^2)/[r])';
'[F]11=[ru]';
'[F]21=[ru]^2/[r]+p';
'[F]31=[ru]*[rv]/[r]';
'[F]41=[ru]*[rw]/[r]';
'[F]51=[ru]*([rE]+p)/[r]';
'[F]12=[rv]';
'[F]22=[ru]*[rv]/[r]';
'[F]32=[rv]^2/[r]+p';
'[F]42=[rv]*[rw]/[r]';
'[F]52=[rv]*([rE]+p)/[r]';
'[F]13=[rw]';
'[F]23=[ru]*[rw]/[r]';
'[F]33=[rv]*[rw]/[r]';
'[F]43=[rw]^2/[r]+p';
'[F]53=[rw]*([rE]+p)/[r]';
};

replacements = {
% U
'[r]',  [U '1'];
'[ru]', [U '2'];
'[rv]', [U '3'];
'[rw]', [U '4'];
'[rE]', [U '5'];
% Fi
'[F]',  Fi;
% GAMMA
'[GAMMA]', GAMMA;
};

for r=replacements'
  out = strrep(out,r{1},r{2});
end


function out = ac_pressure(U, GAMMA, P)
% PRESSURE equation of state
% out = ac_pressure(U, GAMMA, P)
%
% p = (gamma-1)*(rE-1/2*(ru^2+rv^2+rw^2)/r)

if (isnumeric(GAMMA))
  GAMMA = num2str(GAMMA);
end

out = { ...
'[p]=([GAMMA]-1)*([rE]-1/2*([ru]*[ru]+[rv]*[rv]+[rw]*[rw])/[r])'; ...
};

replacements = {
% U
'[r]',  [U '1'];
'[ru]', [U '2'];
'[rv]', [U '3'];
'[rw]', [U '4'];
'[rE]', [U '5'];
% GAMMA
'[GAMMA]', GAMMA;
% P
'[p]', P;
};

for r=replacements'
  out = strrep(out,r{1},r{2});
end


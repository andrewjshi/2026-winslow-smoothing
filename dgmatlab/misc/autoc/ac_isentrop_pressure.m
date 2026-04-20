function out = ac_pressure(U, GAMMA, k_pr, P)
% PRESSURE equation of state
% out = ac_pressure(U, GAMMA, k_pr, P)
%
% p = k_pr*r^GAMMA

if (isnumeric(GAMMA))
  GAMMA = num2str(GAMMA);
end
if (isnumeric(k_pr))
  k_pr = num2str(k_pr);
end

out = { ...
'[p]=[k_pr]*[r]^([GAMMA])'; ...
};

replacements = {
% U
'[r]',  [U '1'];
% GAMMA
'[GAMMA]', GAMMA;
% k_pr
'[k_pr]', k_pr;
% P
'[p]', P;
};

for r=replacements'
  out = strrep(out,r{1},r{2});
end


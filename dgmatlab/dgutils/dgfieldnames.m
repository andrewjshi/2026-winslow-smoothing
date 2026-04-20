function flds=dgfieldnames(type)

switch type
 case {'msh', 'mshsl'}
  flds={'eltype','p','t','t2t','t2n','ncurved','ecurved','tcurved', ...
        'porder','s','s0','tlocal','sbnd','tbndlocal','p1','tpart','partsz','BINs','lt2ts','comm','iBIs'};
 case 'data'
  flds={'eltype','gx','gw','gfs','wgfs','gx0','gw0','gfs0','wgfs0','Nx0','egx','egw','efs','wefs',...
        'egix','permmap','permegix','lM0','lM0R','lM0H','lMref','lMrefR','lMrefH','KWP','KWPsizes','KWPegix'};
 case 'phys'
  flds={'bndcnds','uinf','pars','viscous','sensor_pars','time','ale_pars'};
 otherwise
  error('Unknown structure type.');
end

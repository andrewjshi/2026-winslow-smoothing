function str=mkwrapper(funcname,pars,switchvar,mappings)

str='';
str=[str,sprintf('void %s(',funcname)];
for i=1:numel(switchvar)
  if i==1, pre=''; else pre=', '; end
  str=[str,pre,sprintf('int %s',switchvar{i})];
end

cpars=pars{1};
for j=1:numel(cpars)
  str=[str,sprintf(', double *%s',cpars{j})];
end

str=[str,sprintf(')\n')];
str=[str,sprintf('{\n')];
for i=1:size(mappings,1)
  cparlist='';
  cpars=pars{mappings{i,3}};
  for j=1:numel(cpars)
    if j>=2, cparlist=[cparlist,',']; end
    cparlist=[cparlist,cpars{j}];
  end

  cond='';
  for j=1:numel(switchvar)
    if j==1, pre=''; else pre=' && '; end
    cond=[cond,pre,sprintf('%s==%d',switchvar{j},mappings{i,1}(j))];
  end
  
  if i==1, elsestr=''; else elsestr='else '; end
  
  str=[str,sprintf('  %sif (%s)\n    %s(%s);\n', ...
                   elsestr,cond,mappings{i,2},cparlist)];
end
str=[str,sprintf('}\n')];

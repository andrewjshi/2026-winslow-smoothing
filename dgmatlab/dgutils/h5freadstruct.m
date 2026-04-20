function strct=h5freadstruct(fname,expected_flds)

fid = H5F.open(fname, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

nflds = H5G.get_info(fid).nlinks;
flds = cell(nflds, 1);
for i=1:nflds
  flds{i} = H5G.get_objname_by_idx(fid, i-1);
end

strct=struct;
for ii=1:nflds

  dataset = H5D.open(fid, flds{ii}, 'H5P_DEFAULT');
  space = H5D.get_space(dataset);
  H5S_NULL = 2; %H5ML.get_constant_value('H5S_NULL');
  if H5S.get_simple_extent_type(space) == H5S_NULL
    dims_attr = H5A.open(dataset, 'dims', 'H5P_DEFAULT');
    %dims_attr_space = H5A.get_space(dims_attr);
    %assert(H5S.get_simple_extent_ndims(dims_attr_space) == 1);
    %[~, ndim] = H5S.get_simple_extent_dims(dims_attr_space);
    dims = H5A.read(dims_attr, 'H5ML_DEFAULT');
    H5A.close(dims_attr);

    intcl_attr = H5A.open(dataset, 'intcl', 'H5P_DEFAULT');
    intcl = H5A.read(intcl_attr, 'H5ML_DEFAULT');
    H5A.close(intcl_attr);

    if intcl == 0
      cl = 'int32';
    elseif intcl == 1
      cl = 'double';
    elseif intcl == 2
      cl = 'logical';
    else
      error('Unknown data type');
    end

    if all(dims)
        dat = zeros(dims, cl);
    else
        dat = zeros([0,0], cl);
    end
  else
    dat = H5D.read(dataset, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', ...
                   'H5P_DEFAULT');
    if strcmp(class(dat), 'uint8')
      dat = logical(dat);
    end
  end
  H5S.close(space);
  H5D.close(dataset);

  strct=setfield(strct, flds{ii}, dat);
end

H5F.close(fid);

if nargin >= 2 && 0
  fprintf('Found:');
  fprintf(' %s', flds{:});
  fprintf('\n');
  fprintf('Expected:');
  fprintf(' %s', expected_flds{:});
  fprintf('\n');
end

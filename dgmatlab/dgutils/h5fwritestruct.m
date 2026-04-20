function h5fwritestruct(strct,fname,flds)

%if ~all(ismember(flds,fieldnames(strct)))
%  warning('Not all field names available in structure.');
%end

% These calls to the HDF-5 API should be essentially the same
% as those used in src/iohdf.cpp.
fid = H5F.create(fname, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');

for ii=1:numel(flds)
  if any(ismember(flds{ii},fieldnames(strct)))
    dat=getfield(strct,flds{ii});

    % h5 doesn't know 'logical', so use uint8 instead.
    if islogical(dat)
        dat = uint8(dat);
    end
    ndim = ndims(dat);
    dims = fliplr(size(dat));
    cl = class(dat);
    if strcmp(cl, 'uint8')
        type = 'H5T_NATIVE_UINT8';
        intcl = int32(2);
    elseif strcmp(cl, 'int32')
        type = 'H5T_NATIVE_INT32';
        intcl = int32(0);
    elseif strcmp(cl, 'double')
        type = 'H5T_NATIVE_DOUBLE';
        intcl = int32(1);
    else
        error('Unknown data type');
    end
    dims_type = 'H5T_NATIVE_INT32';
    intcl_type = 'H5T_NATIVE_INT32';

    if numel(dat)
        space = H5S.create_simple(ndim, dims, dims);
        dset = H5D.create(fid, flds{ii}, type, space, 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');
        H5D.write(dset, type, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', dat)
        H5D.close(dset);
        H5S.close(space);
    else
        % Prevent HDF5 error:
        %   "zero sized dimension for non-unlimited dimension"
        %
        % So instead of trying to create an emtpy array, we instead
        % create a NULL dataset with 'dims' and 'intcl' attributes.
        %
        % How stupid.
        space = H5S.create('H5S_NULL');
        dset = H5D.create(fid, flds{ii}, type, space, 'H5P_DEFAULT', 'H5P_DEFAULT', 'H5P_DEFAULT');

        % Write dims
        dims_space = H5S.create_simple(1, [ndim], [ndim]);
        dims_space_attr = H5A.create(dset, 'dims', dims_type, dims_space, 'H5P_DEFAULT');
        H5A.write(dims_space_attr, dims_type, int32(dims));
        H5A.close(dims_space_attr);
        H5S.close(dims_space);

        % Write intcl
        intcl_space = H5S.create('H5S_SCALAR');
        intcl_space_attr = H5A.create(dset, 'intcl', intcl_type, intcl_space, 'H5P_DEFAULT');
        H5A.write(intcl_space_attr, intcl_type, intcl);
        H5A.close(intcl_space_attr);
        H5S.close(intcl_space);

        H5D.close(dset);
        H5S.close(space);
    end
  end
end

H5F.close(fid)

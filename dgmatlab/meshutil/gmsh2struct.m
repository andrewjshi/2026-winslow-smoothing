function gmsh = gmsh2struct(fname)

    fid = fopen(fname,'r');
    if fid == -1
        error('Can''t open file');
    end

    readuntil(fid, '$MeshFormat');
    str = fgetl(fid);
    if strcmp(str, '2.2 0 8') == 0
        fclose(fid);
        error('MeshFormat must be ''2.2 0 8''');
    end
    readuntil(fid, '$EndMeshFormat');

    next_section = fgetl(fid);
    if strcmp(next_section, '$PhysicalNames')
        nbrnames = fscanf(fid, '%d', 1);
        names = cell(nbrnames,3);
        for ii = 1:nbrnames
            names{ii,1} = fscanf(fid, '%d', 1);
            names{ii,2} = fscanf(fid, '%d ', 1);
            names{ii,3} = fgetl(fid);
        end
        gmsh.PhysicalNames = names;
        readuntil(fid, '$EndPhysicalNames');
        next_section = fgetl(fid);
    else
        gmsh.PhysicalNames = {};
        if ~strcmp(next_section, '$Nodes')
            error('No $Nodes section');
        end
    end

    np = fscanf(fid, '%d', 1);
    p = zeros(np, 3);
    for ii = 1:np
        nbr = fscanf(fid, '%d', 1);
        if (nbr ~= ii)
            fclose(fid);
            error('Only supports contiguous node numbering');
        end
        p(ii,:) = fscanf(fid, '%f', 3)';
        fgetl(fid);
    end
    gmsh.Nodes = p;
    
    readuntil(fid, '$Elements');
    ne = fscanf(fid, '%d', 1);
    e.type = zeros(ne,1);
    e.tags = zeros(ne,2);
    nn = cell(ne,1);
    max_nn = 0;
    for ii = 1:ne
        nbr = fscanf(fid, '%d', 1);
        e.type(ii) = fscanf(fid, '%d', 1);
        nbrtags = fscanf(fid, '%d', 1);
        if (nbrtags ~= 2)
            fclose(fid);
            error('Only supports exactly 2 tags');
        end
        e.tags(ii,:) = fscanf(fid, '%d', nbrtags);
        node_number_str = fgetl(fid);
        nn{ii} = sscanf(node_number_str, '%d', inf);
        max_nn = max(max_nn, numel(nn{ii}));
    end
    e.node_numbers = zeros(ne, max_nn);
    for ii = 1:ne
        e.node_numbers(ii,1:numel(nn{ii})) = nn{ii};
    end
    gmsh.Elements = e;
    
    fclose(fid);
end

function readuntil(fid,str)
    found = false;
    while ~feof(fid)
        fline = fgetl(fid);
        if ~isempty(fline) & ~isempty(strmatch(str, fline))
            found = true;
            break;
        end
    end
    if ~found
        fclose(fid);
        error(['File error - cannot find string ', str]);
    end
end

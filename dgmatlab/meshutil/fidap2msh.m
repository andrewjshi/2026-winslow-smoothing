function msh = fidap2msh(fname)

    fid = fopen(fname,'r');
    if fid == -1
        error('Can''t open file');
    end

    readuntil(fid, 'FIDAP NEUTRAL FILE');
    str = fgetl(fid);
    str = fgetl(fid);
    str = fgetl(fid);
    str = fgetl(fid);

    nn = fscanf(fid, '%d', 2);
    np = nn(1);
    nel_total = nn(2);
    
    readuntil(fid, 'NODAL COORDINATES');
    
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

    readuntil(fid, 'ELEMENT GROUPS');
    
    cgroup = 0;
    while ~feof(fid)
        cgroup = cgroup + 1;
        nn = fscanf(fid, 'GROUP: %d ELEMENTS: %d NODES: %d', 3);
        fgetl(fid);
        str = fscanf(fid, 'ENTITY NAME: %s', 1);
        group = nn(1);
        nel = nn(2);
        nelnodes = nn(3);
        fgetl(fid);
        
        if group == 1 & nelnodes ~= 4
            error('Only supports group == 1 with tets');
        elseif group > 1 & nelnodes ~= 3
            warning('Non-triangles in group > 1 -- ignoring');
            for iel = 1:nel
                fgetl(fid);
            end
            cgroup = cgroup - 1;
            continue
        end
        
        el{cgroup} = zeros(nel, nelnodes);
        for iel = 1:nel
            nbr = fscanf(fid, '%d', 1);
            el{cgroup}(iel,:) = fscanf(fid, '%d', nelnodes)';
            fgetl(fid);
        end
        groupnames{cgroup} = str;
    end
    t = el{1};
    el = el(2:end);
    bndnames = groupnames(2:end);
    
    msh = ml2msh(p, t);
    
    % Mark boundary tris to make sure we label all of them
    msh.t2t(msh.t2t == -1) = -999; 
    t2t = msh.t2t';
    
    % All faces
    faces = mkfaces(t);
    faces = sort(faces,2);

    % Identify boundary triangles
    for i = 1:numel(el)
        [mem, loc] = ismember(sort(el{i},2), faces, 'rows');
        if ~all(mem)
            error('Boundary triangle not found in mesh');
        end
        
        if ~all(t2t(loc)) == -999
            error('Boundary triangle not an external face in mesh');
        end
        
        t2t(loc) = -i;
    end

    if any(t2t(:) == -999)
        error('Some external faces in mesh not marked as boundary groups');
    end

    msh.t2t = t2t';
    msh.bndnames = bndnames;
    
    fclose(fid);
end

function readuntil(fid,str)
    found = false;
    while ~feof(fid)
        fline = fgetl(fid);
        if ~isempty(fline) & ~isempty(findstr(str, fline))
            found = true;
            break;
        end
    end
    if ~found
        fclose(fid);
        error(['File error - cannot find string ', str]);
    end
end

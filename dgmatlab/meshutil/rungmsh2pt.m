function [p,t] = rungmsh2pt(gmshfname, varargin)
% Run gmsh on input filename, read output and return p,t mesh.

    fout = [tempname, '.msh'];
    gmshstr = ['gmsh -3 -format msh2 -o ', fout, ' ', gmshfname];
    for arg = varargin
        gmshstr = [gmshstr, ' ', arg{1}];
    end
    [status, out] = system(gmshstr);
    if status ~= 0
        error(sprintf('Something went wrong running gmsh: See output below\n\n%s', out));
    end

    gmsh = gmsh2struct(fout);
    delete(fout);
    
    p = gmsh.Nodes;
    
    ntq = [sum(gmsh.Elements.type == 2), sum(gmsh.Elements.type == 3)];
    nth = [sum(gmsh.Elements.type == 4), sum(gmsh.Elements.type == 5)];

    if any(nth > 0) % 3D
        if nth(1) > 0 & nth(2) == 0  % Tets
            t = gmsh.Elements.node_numbers(gmsh.Elements.type == 4, 1:4);
        elseif nth(1) == 0 & nth(2) > 0 % Hexes
            t = gmsh.Elements.node_numbers(gmsh.Elements.type == 5, 1:8);
        else
            error('Cannot handle mixed element 3D meshes.');
        end
    elseif any(ntq > 0) % 2D
        if ntq(1) > 0 & ntq(2) == 0  % Tris
            t = gmsh.Elements.node_numbers(gmsh.Elements.type == 2, 1:3);
        elseif ntq(1) == 0 & ntq(2) > 0 % Quads
            t = gmsh.Elements.node_numbers(gmsh.Elements.type == 3, [1,2,4,3]);
        else
            error('Cannot handle mixed triangle / quad meshes.');
        end
        
        p = p(:,1:2); % Remove z-coordinate since 2D
        
    else
        error('Must have nonzero # of either tris, quads, tets, or hexes.');
    end
end

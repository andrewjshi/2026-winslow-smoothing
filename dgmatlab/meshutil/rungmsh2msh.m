function msh = rungmsh2msh(gmshfname, varargin)
% Run gmsh on input filename, read output and return 3DG mesh.

    opts = struct;
    fout = [tempname, '.msh'];
    gmshstr = ['gmsh -3 -format msh2 -o ', fout, ' ', gmshfname];
    for arg = varargin
        if ischar(arg{1})
            gmshstr = [gmshstr, ' ', arg{1}];
        elseif isstruct(arg{1})
            opts = arg{1};
        else
            error('Syntax error.');
        end
    end
    [status, out] = system(gmshstr);
    if status ~= 0
        error(sprintf('Something went wrong running gmsh: See output below\n\n%s', out));
    end

    msh = gmsh2msh(fout, opts);
    delete(fout);
end

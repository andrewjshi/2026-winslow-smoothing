function save_gmsh_points(filename, pv)
    fid = fopen(filename, 'w');
    for i = 1:size(pv,1)
        fprintf(fid, 'Point(%d) = {%g, %g, 0};\n', i, pv(i,:));
    end
    fclose(fid);
end

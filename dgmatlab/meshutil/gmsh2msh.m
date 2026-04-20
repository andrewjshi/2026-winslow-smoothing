function [msh,gmsh] = gmsh2msh(gmsh_fname, opts)
%GMSH2MSH  Load GMSH mesh file and convert to 3DG msh format
%
%   [MSH,GMSH] = GMSH2MSH(GMSH_FNAME)
    
    if nargin<2, opts = struct; end
    if ~isfield(opts,'intbnd'), opts.intbnd = false; end
    
    %%%%%%%%%%%%%%
    % GMSH data
    
    % Element types
    lins = [1,8,26,27,28];
    tris = [2,9,21,23,25];
    quads = [3,10,36,37,38];
    tets = [4,11,29,30,31];
    hexs = [5,12,92,93,94];
    % Dimensions
    eltype2dim = zeros(1,100);
    eltype2dim(lins) = 1;
    eltype2dim(tris) = 2;
    eltype2dim(quads) = 2;
    eltype2dim(tets) = 3;
    eltype2dim(hexs) = 3;
    
    % Node ordering
    tri_node_order_map = {[3,1,2];
                          [3,6,1,5,4,2];
                          [3,8,9,1,7,10,4,6,5,2];
                          [3,10,11,12,1,9,15,13,4,8,14,5,7,6,2];
                          [3;12;13;14;15;1;11;18;21;16;4;10;20;19;5;9;17;6;8;7;2]
                        };

    tet_node_order_map = {[4,1,2,3];
                          [4,8,1,10,5,2,9,7,6,3];
                          [4;11;12;1;15;18;5;16;6;2;13;19;10;20;17;7;14;9;8;3];
                          [4;14;15;16;1;20;28;26;5;21;27;6;22;7;2;17;30;29;13;32;35;23;33;25;8;18;31;12;34;24;9;19;11;10;3];
                          [4;17;18;19;20;1;25;37;40;35;5;26;39;38;6;27;36;7;28;8;2;21;42;44;41;16;47;56;53;29;50;54;34;48;31;9;22;45;46;15;52;55;32;51;33;10;23;43;14;49;30;11;24;13;12;3]
                        };    

    quad_node_order_map = {[1,2,4,3];
                           [1,5,2,8,9,6,4,7,3];
                           [1;5;6;2;12;13;14;7;11;16;15;8;4;10;9;3];
                           [1;5;6;7;2;16;17;21;18;8;15;24;25;22;9;14;20;23;19;10;4;13;12;11;3];
                           [1;5;6;7;8;2;20;21;25;26;22;9;19;32;33;34;27;10;18;31;36;35;28;11;17;24;30;29;23;12;4;16;15;14;13;3]
                        };
    
    hex_node_order_map = {[1,2,4,3,5,6,8,7];
                          [1;9;2;10;21;12;4;14;3;11;22;13;23;27;24;16;25;15;5;17;6;18;26;19;8;20;7];
                          [1;9;10;2;11;33;36;15;12;34;35;16;4;20;19;3;13;37;38;17;41;57;58;45;44;60;59;46;23;50;49;21;14;40;39;18;42;61;62;48;43;64;63;47;24;51;52;22;5;25;26;6;27;53;54;29;28;56;55;30;8;32;31;7];
                          [1;9;10;11;2;12;45;52;48;18;13;49;53;51;19;14;46;50;47;20;4;26;25;24;3;15;54;58;55;21;63;99;107;100;72;70;108;119;110;76;66;102;112;101;73;30;82;85;81;27;16;61;62;59;22;67;109;120;111;79;71;121;125;122;80;69;114;123;113;77;31;86;89;88;28;17;57;60;56;23;64;103;115;104;75;68;116;124;117;78;65;106;118;105;74;32;83;87;84;29;5;33;34;35;6;36;90;94;91;39;37;97;98;95;40;38;93;96;92;41;8;44;43;42;7];
                          [1;9;10;11;12;2;13;57;68;67;60;21;14;61;69;72;66;22;15;62;70;71;65;23;16;58;63;64;59;24;4;32;31;30;29;3;17;73;77;78;74;25;89;153;161;162;154;105;100;163;185;188;167;109;99;164;186;187;168;110;92;156;172;171;155;106;37;122;126;125;121;33;18;84;85;86;79;26;93;165;189;190;169;116;101;193;209;210;197;117;104;196;212;211;198;118;98;175;202;201;173;111;38;127;134;133;132;34;19;83;88;87;80;27;94;166;192;191;170;115;102;194;213;214;200;120;103;195;216;215;199;119;97;176;203;204;174;112;39;128;135;136;131;35;20;76;82;81;75;28;90;157;177;178;158;108;95;179;205;206;181;114;96;180;208;207;182;113;91;160;184;183;159;107;40;123;129;130;124;36;5;41;42;43;44;6;45;137;141;142;138;49;46;148;149;150;143;50;47;147;152;151;144;51;48;140;146;145;139;52;8;56;55;54;53;7]};
    
    %%%%%%%%%%%%%%

    gmsh = gmsh2struct(gmsh_fname);
    
    dims = unique(eltype2dim(gmsh.Elements.type));
    dim = max(dims);
    eltypes = unique(gmsh.Elements.type);
    p = gmsh.Nodes(:,1:dim);
    
    switch dim
      case 2
        is_tri = any(ismember(tris, eltypes));
        is_quad = any(ismember(quads, eltypes));
        if is_tri & is_quad
            error('Cannot handle mixed tri/quad meshes');
        elseif is_tri
            vols = tris;
            surfs = lins;
            node_order_map = tri_node_order_map;
            tmap = 1:3;
            eltype = t_simplex;
        elseif is_quad
            vols = quads;
            surfs = lins;
            node_order_map = quad_node_order_map;
            tmap = [1,2,4,3];
            eltype = t_block;
        else
            error('Must have either tris or quads in 2D');
        end
      case 3
        is_tet = any(ismember(tets, eltypes));
        is_hex = any(ismember(hexs, eltypes));
        if is_tet & is_hex
            error('Cannot handle mixed tet/hex meshes');
        elseif is_tet
            vols = tets;
            surfs = tris;
            node_order_map = tet_node_order_map;
            tmap = 1:4;
            eltype = t_simplex;
        elseif is_hex
            vols = hexs;
            surfs = quads;
            node_order_map = hex_node_order_map;
            tmap = [1,2,4,3,5,6,8,7];
            eltype = t_block;
        else
            error('Must have either tets or hexes in 3D');
        end
      otherwise
        error('Only implemented in 2D and 3D');
    end

    porder = find(ismember(vols, eltypes));
    if length(porder) ~= 1
        error('Must have unique order volume elements');
    end
    
    fprintf('  Dimension %d, porder %d\n\n', dim, porder);

    t1 = gmsh.Elements.node_numbers(gmsh.Elements.type == vols(porder), :);
    t = t1(:, tmap);
    [p_straight, t_straight] = fixmesh(p,t);
    msh = ml2msh(p_straight, t_straight);
    msh = mshcurved(msh, 'all');
    
    if isfield(gmsh, 'PhysicalNames')
        names = gmsh.PhysicalNames;
        if ~isempty(names)
            names = names(cell2mat(names(:,1)) == dim - 1, 2:3);
            bndnbrs = cell2mat(names(:,1));
            names = names(:,2);
            msh.bndnames = cell(1, max(bndnbrs));
            for ii = 1:numel(bndnbrs)
                msh.bndnames{bndnbrs(ii)} = names{ii};
            end
        end
    end

    if any(gmsh.Elements.tags(:,1) == 0) % Use elementary entities
        tagcol = 2;
    else                                 % Use physical entities
        tagcol = 1;
    end
    
    surf_porder = find(ismember(surfs, eltypes));
    if  ~isempty(surf_porder)
        if length(surf_porder) > 1
            error('Must have unique order surface elements');
        end
        if porder ~= surf_porder
            error('Must have same porder for volume and surface elements');
        end
        
        % Assign boundary numbers
        faces = mkfaces(t, eltype);
        surf_loc = gmsh.Elements.type == surfs(porder);
        surf = gmsh.Elements.node_numbers(surf_loc, 1:size(faces,2));
        nbr = gmsh.Elements.tags(surf_loc, tagcol);

        [member, ix] = ismember(sort(surf,2), sort(faces,2), 'rows');
        if ~all(member)
            %warning('Surface elements not found in volume elements');
            ix = ix(member);
            nbr = nbr(member);
        end
        t2t = msh.t2t';
        if ~all(t2t(ix) == -1) & ~opts.intbnd
            % Remove interior boundaries
            ixx = t2t(ix) == -1;
            ix = ix(ixx);
            nbr = nbr(ixx);
        end
        t2t(ix) = -nbr;
        msh.t2t = t2t';
    end
    
    % High-order nodes
    msh = nodealloc(msh, (0:porder)'/porder);
    t1map = t1(:,node_order_map{porder});
    msh.p1 = permute(reshape(p(t1map(:), :), [size(t1map), dim]), [2,3,1]);

end

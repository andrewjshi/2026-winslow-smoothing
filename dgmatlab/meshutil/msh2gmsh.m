function msh2gmsh(fname, bndnames, msh)
%MSH2GMSH
%   Syntax:
%      msh2gmsh(fname, bndnames);
%      msh2gmsh(fname, bndnames, msh);
    
    [fpath,fname,fext] = fileparts(fname);
    fname_out = [fullfile(fpath, [fname, '.msh'])];
    
    if nargin < 3
        msh = h5freadstruct(fname);
    end

    nbndtypes = double(-min(msh.t2t(:)));

    if nargin<2 | isempty(bndnames)
        bndnames = num2cell(1:nbndtypes);
    else
        if nbndtypes ~= numel(bndnames)
            error('Number of bnd names does not match with msh.t2t');
        end
    end
    
    [ns,dim,nt] = size(msh.p1);
    data = dginit(msh);
    porder = double(msh.porder);
    eltype = msh.eltype;
    
    %%%%%%%%%%%%%%
    % GMSH data
    
    % Element types
    lins = [1,8,26,27,28];
    tris = [2,9,21,23,25];
    quads = [3,10,36,37,38];
    tets = [4,11,29,30,31];
    hexs = [5,12,92,93,94];

    % Node ordering
    lin_node_order_map = {[1,2];
                          [1,3,2];
                          [1,4,2,3];
                          [1,5,2,3,4];
                          [1,6,2,3,4,5]
                        };
    
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

    p1 = reshape(permute(msh.p1, [1,3,2]), [], dim);
    snap = eps * 2^20;
    [~,ix,jx] = unique(round(p1/snap) * snap, 'rows');
    p1 = p1(ix, :);
    np1 = size(p1,1);

    if dim == 2
        p1(:,3) = 0.0;
    end
    
    elems = reshape(jx, ns, nt);

    % gmsh_node_map = {{[0,1], [0,2,1], [0,2,3,1], [0,2,3,4,1]},
    %                  {[0,1,2], [0,3,1,5,4,2], [0,3,4,1,8,9,5,7,6,2], ...
    %                   [0,3,4,5,1,11,12,13,6,10,14,7,9,7,2]}};
    % gmsh_node_imap = cell(size(gmsh_node_map));
    % for i = 1:numel(gmsh_node_map)
    %     for j = 1:numel(gmsh_node_map{i})
    %         [~,gmsh_node_imap{i}{j}] = sort(gmsh_node_map{i}{j});
    %     end
    % end

    switch char([dim, eltype])
      case char([2, t_simplex])
        gmsh_element_type = { lins, tris};
        gmsh_node_imap = { lin_node_order_map, tri_node_order_map };
      case char([2, t_block])
        gmsh_element_type = { lins, quads };
        gmsh_node_imap = { lin_node_order_map, quad_node_order_map };
      case char([3, t_simplex])
        gmsh_element_type = { tris, tets};
        gmsh_node_imap = { tri_node_order_map, tet_node_order_map };
      case char([3, t_block])
        gmsh_element_type = { quads, hexs};
        gmsh_node_imap = { quad_node_order_map, hex_node_order_map };
      otherwise
        error('Unsupported elements');
    end

    
    nbndelems = 0;
    bndelems = cell(1, nbndtypes);
    for it = 1:nt
        for j = 1:size(data.egix,2)
            jt = msh.t2t(j,it);
            if jt < 0
                cbnd = -jt;
                cbndelem = elems(data.egix(:,j)+1, it);
                bndelems{cbnd}(end+1,:) = cbndelem;
                nbndelems = nbndelems + 1;
            end
        end
    end

f = fopen(fname_out, 'w');

fprintf(f, '$MeshFormat\n');
fprintf(f, '2.2 0 8\n');
fprintf(f, '$EndMeshFormat\n');

fprintf(f, '$PhysicalNames\n');
fprintf(f, '%d\n', nbndtypes + 1);
for i = 1:nbndtypes
    fprintf(f, '%d %d "%s"\n', dim-1, i, bndnames{i});
end
fprintf(f, '%d %d "Interior"\n', dim, nbndtypes + 1);
fprintf(f, '$EndPhysicalNames\n');

fprintf(f, '$Nodes\n');
fprintf(f, '%d\n', np1);
fprintf(f, '%d %.16g %.16g %.16g\n', [1:np1; p1']);
fprintf(f, '$EndNodes\n');

fprintf(f, '$Elements\n');
fprintf(f, '%d\n', nt + nbndelems);

cpos = 0;
[~,map] = sort(gmsh_node_imap{1}{porder});
gmsh_eltype = gmsh_element_type{1}(porder);
for i = 1:numel(bndelems)
    cbnd = bndelems{i};
    [ne,nse] = size(cbnd);
    if ne > 0
        o = ones(1,ne);
        to_gmsh = [(1:ne) + cpos; gmsh_eltype*o; 2*o; i*[o;o]; cbnd(:,map)'];
        fprintf(f, [repmat('%d ', 1, nse+5), '\n'], to_gmsh);
        cpos = cpos + ne;
    end
end

o = ones(1,nt);
[~,map] = sort(gmsh_node_imap{2}{porder});
gmsh_eltype = gmsh_element_type{2}(porder);
to_gmsh = [cpos+(1:nt); gmsh_eltype*o; 2*o; (nbndtypes + 1)*[o;o]; elems(map,:)];
fprintf(f, [repmat('%d ', 1, ns+5), '\n'], to_gmsh);

fprintf(f, '$EndElements\n');

fclose(f);

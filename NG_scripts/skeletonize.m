%% skeletonize
% Skeletonization function using TEASAR algorithm to build the skeleton of the cell

function [nodes, edges, node_radii, rt_coord, dt_coord] = skeletonize(segmentation, object_id, scale, const, source)
    %% Inputs
    % im_points - 3D array segmentation chunk
    % input_resolution - self explanatory, difference in real length between
    %                    point [1 1 1] and [0 0 0]
    % scale and const - "scale" and "const" parameter from TEASAR paper,
    %                   larger values mean fewer hairs
    % source - n x 3 array of source coordinates
    
    %% Ouputs
    % nodes - n x 3 array of node coordintes
    % edges - n x 2 array of node numbers (col1: start, col2: end) 
    % root_node - Node number for the root node
    % node_radii - DBF at the nodes
    % rt_coord - n x 3 array of root coordinates
    % dt_coord - n x 3 array of destination coordinates
    
    obj_points = array2ind(segmentation,object_id);
    
    % Initial inputs
  if ~exist('object_id','var') || isempty(object_id)
        object_id = [];
    end
    if ~exist('scale','var') || isempty(scale)
        scale = 6;
    end
    if ~exist('const','var') || isempty(const)
        const = 6;
    end
    if ~exist('source','var') || isempty(source)
        source = obj_points(1,:);
    end
    
    
    tic
    disp('Building skeleton...');
    [cube_out, cube_paths] = cube_TEASAR(obj_points, scale, const, source);
    toc
    
    rt = cube_out.root_node;
    rt_coord = cube_out.nodes(rt,:);
    dt = cube_out.dest;
    dt_coord = cube_out.nodes(dt,:);
  
    disp('Source Node');
    disp(rt_coord); % Display root node coordinates
    disp('Destination Node'); 
    disp(dt_coord); % Display destination node coordinates
                      
        
    % Consolidate nodes and edges
    num_nodes = length(cube_out.node_radii);
    num_edges = size(cube_out.edges,1);
    
    nodes = zeros(num_nodes,3);
    node_radii = zeros(num_nodes,1);
    edges = zeros(num_edges,2);
        

    n = single(cube_out.nodes);
    nodes(1:num_nodes, :) = n;
    node_radii(1:num_nodes) = cube_out.node_radii;
    edges(1:num_edges, :) = cube_out.edges;
    
    
    disp('resolving edges'); tic;
    dup_nodes = find_duplicate_entries(nodes);    
    for n = size(dup_nodes,1):-1:1;
        edges(edges==dup_nodes(n,2)) = dup_nodes(n,1);
    end
    node_list = unique(edges(:));
    
    nodes = nodes(node_list,:);
    node_radii = node_radii(node_list);
    for n = 1:length(node_list)
        edges(edges==node_list(n)) = n;
    end
    edges = unique(edges, 'rows');
    toc
       
end


function [cube_out, cube_paths] = cube_TEASAR(p, scale, const, source)
    NDIMS = size(p,2);

    quick_ind = @(S, I) (sub2ind(S, I(:,1), I(:,2), I(:,3)));
    im_max = max(p);
    
    bin_im = false(im_max);
    
    bin_im(quick_ind(im_max, p)) = true;
    node2ind = find(bin_im(:));
    loc2node = zeros(size(bin_im));
    loc2node(node2ind) = 1:length(node2ind);
    
    N = length(node2ind);
    disp('Number of points')
    disp(N);
    
    ind2node = sparse(node2ind,ones(N,1),(1:N)',prod(im_max),1);
    
    node2loc = zeros(N, 3);
    [node2loc(:,1), node2loc(:,2), node2loc(:,3)] = ind2sub(im_max, node2ind); 

    
    % Distance to the boundary
    DBF = bwdist(~bin_im); 
    
    
    % Penalty weight for the edges 
    M = max(DBF(:))*1.01;
    p_v = 5000*(1-DBF/M).^16;

    

    % Create graph, use euclidean distance for edge weights
    % 26-connectivity
    [nhood(:,1) nhood(:,2) nhood(:,3)] = ind2sub([3 3 3], 1:27);
    nhood = nhood - 2;
    nhood(all(nhood==0,2),:) = [];
    hood_weight = sqrt(sum(nhood.^2,2));
    hood_length = size(nhood,1);
    
    dest_node = zeros(N,hood_length);
    my_node = zeros(N,hood_length);
    edge_weight = zeros(N,hood_length);
    edge_dist = zeros(N,hood_length);
    
    
    for n = 1:hood_length
        my_node(:,n) = 1:N;        
        
        dest_loc = node2loc + ones(N,1)*nhood(n,:);

        is_valid = all(dest_loc >= ones(N,1)*min(p), 2) & all(dest_loc <= ones(N,1)*im_max,2);
        dest_ind = zeros(N,1);
        dest_ind(is_valid) = sub2ind(size(bin_im), dest_loc(is_valid,1), dest_loc(is_valid,2), dest_loc(is_valid,3));
        dest_node(is_valid, n) = ind2node(dest_ind(is_valid));
        is_valid = dest_node(:,n) ~= 0;
        edge_dist(is_valid, n) = hood_weight(n);
        edge_weight(is_valid, n) = hood_weight(n) * p_v(node2ind(dest_node(is_valid,n))); 
    end

  
    has_edge = dest_node(:) ~= 0;
    G = sparse(my_node(has_edge), dest_node(has_edge), edge_weight(has_edge), N, N); % Graph with euclidean distance
    G_dist = sparse(my_node(has_edge), dest_node(has_edge), edge_dist(has_edge), N, N); % Graph with penalty
    

    if ~exist('source','var') || isempty(source)
        root_ind = find(bin_im(:),1,'first');
    else
        root_ind = quick_ind(im_max, source);
    end
    

    root_nodes = ind2node(root_ind);
    
      
    cube_out.nodes = [];
    cube_out.edges = [];
    cube_out.root_node = [];
    cube_out.node_radii = [];
    cube_out.dest = [];
    
    total_path = 1;
    is_disconnected = true(1, size(G,1));
    n_root = length(root_nodes);
    t = 1;
    while any(is_disconnected)
        if t <= n_root
            root_node = root_nodes(t);
        
            PDRF_G = graphshortestpath(G_dist,root_node);
            is_disconnected = is_disconnected & isinf(PDRF_G);
            PDRF = zeros(im_max);
            PDRF(node2ind) = PDRF_G;
        else
            root_node = find(is_disconnected,1,'first');
            
            PDRF_G = graphshortestpath(G_dist,root_node);
            is_disconnected = is_disconnected & isinf(PDRF_G);
            PDRF = zeros(im_max);
            PDRF(node2ind) = PDRF_G;
        end
        
        paths = [];
        k = 1;
        mask_im = bin_im & ~isinf(PDRF);
        
        [~, part_paths] = graphshortestpath(G,root_node);
        while any(mask_im(:))
            mask_inds = find(mask_im(:));

            [~, farthest_in_ind] = max(PDRF(mask_inds));
            farthest_node = ind2node(mask_inds(farthest_in_ind));
            paths{k} = part_paths{farthest_node};
%             disp('Path added');
            
            
            for c = 1:length(paths{k})
                my_ind = node2ind(paths{k}(c));
                [my_loc(1) my_loc(2) my_loc(3)] = ind2sub(im_max, my_ind);

                r = DBF(my_ind)*scale + const;
                
                r_begin = ceil(my_loc - r);
                r_begin(r_begin<1) = 1;

                r_end = floor(my_loc + r);
                r_end = min([r_end; im_max]);

                mask_im(r_begin(1):r_end(1), r_begin(2):r_end(2), r_begin(3):r_end(3)) = false;
            end
            k = k + 1;
        end          
        t = t + 1;
   

        % Consolidate paths
        tree_nodes = [];
        for path = 1:length(paths);
            tree_nodes = [tree_nodes; paths{path}'];
        end
        tree_nodes = unique(tree_nodes);
        tree_size = length(tree_nodes);

        inv_tree_nodes = sparse(tree_nodes, ones(tree_size,1), (1:tree_size)', N,1);
        root_node = inv_tree_nodes(paths{1}(1)) + size(cube_out.nodes,1);

        edges = []; 
        dest = [];
        for path = 1:length(paths);
            edges = [edges; full([inv_tree_nodes(paths{path}(1:end-1)') inv_tree_nodes(paths{path}(2:end)')])];
            dest = [dest; inv_tree_nodes(paths{path}(end))+size(cube_out.nodes,1)];
            cube_paths{total_path} = inv_tree_nodes(paths{path}(1:end)');
            cube_paths{total_path} = cube_paths{total_path} + size(cube_out.nodes,1);
            total_path = total_path + 1;
        end
        edges = unique(edges, 'rows');


        nodes = [];
        [nodes(:,1) nodes(:,2) nodes(:,3)] = ind2sub(im_max, node2ind(tree_nodes));
        node_radii = zeros(size(nodes,1),1, 'single');
        for n = 1:size(nodes,1)
            node_radii(n) = DBF(nodes(n,1), nodes(n,2), nodes(n,3));
        end
        
        cube_out.edges = [cube_out.edges; edges + size(cube_out.nodes,1)];
        cube_out.nodes = [cube_out.nodes; nodes];
        cube_out.root_node = [cube_out.root_node; root_node];
        cube_out.node_radii = [cube_out.node_radii; node_radii];
        cube_out.dest = [cube_out.dest; dest];
        
    end
end


function dup_nodes = find_duplicate_entries(nodes)    
    
    [sorted_nodes, sort_inds] = sortrows(nodes);
    
    dup_nodes = zeros(size(nodes,1), 2);
    k = 1;
    for n = 1:size(sorted_nodes,1)-1
        if all(sorted_nodes(n,:) == sorted_nodes(n+1,:))
            dup_nodes(k,1:2) = [n, n+1];
            k = k+1;
        end
    end
    dup_nodes = dup_nodes(1:k-1,:);
    dup_nodes = sort_inds(dup_nodes);
    if size(dup_nodes,2) == 1
        dup_nodes = dup_nodes';
    end
    
    dup_nodes(dup_nodes(:,1) > dup_nodes(:,2),:) = dup_nodes(dup_nodes(:,1) > dup_nodes(:,2),[2 1]);
    
    dup_nodes = sortrows(dup_nodes);
end


function point_cloud = array2ind(array,segm_id)

    if ~exist('segm_id','var') || isempty(segm_id)
        [point_cloud(:,1), point_cloud(:,2), point_cloud(:,3)] = ind2sub(size(array),find(array));
    else
        [point_cloud(:,1), point_cloud(:,2), point_cloud(:,3)] = ind2sub(size(array),find(array==segm_id));
    end
    
end
    
      
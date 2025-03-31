clc; clear; close all;

%% Step 1: Draw an Arbitrary Shape
disp('Draw an arbitrary closed shape by clicking points. Press Enter to finish.');
figure;
axis([0 10 0 10]); % Adjust axis limits as needed
hold on; grid on;
xlabel('X'); ylabel('Y');
title('Draw an Arbitrary Shape (Click Points, Press Enter to Finish)');

[x_shape, y_shape] = ginput;  % User clicks to define the shape
x_shape = [x_shape; x_shape(1)];  
y_shape = [y_shape; y_shape(1)];
plot(x_shape, y_shape, 'b-', 'LineWidth', 2);
scatter(x_shape, y_shape, 50, 'r', 'filled');

%% Step 2: Find and Extend Bounding Box
min_x = min(x_shape);
max_x = max(x_shape);
min_y = min(y_shape);
max_y = max(y_shape);

box_width = max_x - min_x;
min_x_extended = min_x - box_width;
max_x_extended = max_x + box_width;

rectangle('Position', [min_x_extended, min_y, max_x_extended - min_x_extended, max_y - min_y], 'EdgeColor', 'k', 'LineWidth', 2);
disp(['Extended Bounding Box: Min X = ', num2str(min_x_extended), ', Max X = ', num2str(max_x_extended)]);

%% Step 3: User Input for Mesh and Angle
num_points = input('Enter the number of equidistant points on the upper boundary: ');
theta = input('Enter the angle (in degrees) from the vertical: ');
theta_rad = deg2rad(theta);  

upper_boundary_x = linspace(min_x_extended, max_x_extended, num_points);
upper_boundary_y = max_y * ones(1, num_points);
scatter(upper_boundary_x, upper_boundary_y, 50, 'r', 'filled');

% Label upper boundary points with odd numbers
upper_labels = 1:2:2*num_points;
for i = 1:num_points
    text(upper_boundary_x(i), upper_boundary_y(i)+0.2, num2str(upper_labels(i)), 'FontSize', 10, 'FontWeight', 'bold');
end

%% Step 4: User Input for Bridge Length
bridge_length = input('Enter the bridge length (extrusion distance) in mm: ');

%% Step 5: Generate 2D Mesh with Intersection Points
lines = [];
line_from_points = []; % Store which boundary point each line came from

for i = 1:num_points
    x_start = upper_boundary_x(i);
    y_start = upper_boundary_y(i);

    x_end_acute = x_start + (max_y - min_y) * tan(theta_rad);  
    y_end_acute = min_y;  
    x_end_obtuse = x_start - (max_y - min_y) * tan(theta_rad);
    y_end_obtuse = min_y;  

    if x_end_acute > max_x_extended
        x_end_acute = max_x_extended;
        y_end_acute = y_start - (x_end_acute - x_start) / tan(theta_rad);
    end
    if x_end_obtuse < min_x_extended
        x_end_obtuse = min_x_extended;
        y_end_obtuse = y_start - (x_start - x_end_obtuse) / tan(theta_rad);
    end

    m_acute = (y_end_acute - y_start) / (x_end_acute - x_start);
    c_acute = y_start - m_acute * x_start; 
    lines = [lines; m_acute, c_acute, x_start, y_start, x_end_acute, y_end_acute];
    line_from_points = [line_from_points; i, 0]; % 0 for acute angle line 

    m_obtuse = (y_end_obtuse - y_start) / (x_end_obtuse - x_start);
    c_obtuse = y_start - m_obtuse * x_start;
    lines = [lines; m_obtuse, c_obtuse, x_start, y_start, x_end_obtuse, y_end_obtuse];
    line_from_points = [line_from_points; i, 1]; % 1 for obtuse angle line

    plot([x_start, x_end_acute], [y_start, y_end_acute], 'g-', 'LineWidth', 1.5); 
    plot([x_start, x_end_obtuse], [y_start, y_end_obtuse], 'm-', 'LineWidth', 1.5); 
end

%% Step 6: Compute Intersection Points for All Lines
disp('Computing all intersection points within the bounding box:');

% Initialize arrays to store all points
all_points = [];
all_intersection_info = []; % Store boundary point number and intersection info
current_label = 2; % Start with even number 2 for intersection points

% Process in order of upper boundary points with odd numbers (1, 3, 5, etc.)
for boundary_idx = 1:num_points
    % Current boundary point label
    boundary_label = upper_labels(boundary_idx);
    disp(['Processing boundary point ', num2str(boundary_label), ':']);
    
    % Add the boundary point to all_points
    x_boundary = upper_boundary_x(boundary_idx);
    y_boundary = upper_boundary_y(boundary_idx);
    all_points = [all_points; x_boundary, y_boundary];
    all_intersection_info = [all_intersection_info; boundary_idx, 0, 0, boundary_label]; % 0 indicates boundary point
    disp(['  Boundary point ' num2str(boundary_label) ': (', num2str(x_boundary), ', ', num2str(y_boundary), ')']);
    
    % Find the acute angle line from this boundary point (green line)
    acute_line_idx = find(line_from_points(:,1) == boundary_idx & line_from_points(:,2) == 0);
    
    if ~isempty(acute_line_idx)
        % Get line parameters
        line_i = lines(acute_line_idx,:);
        m1 = line_i(1);
        c1 = line_i(2);
        x_start_i = line_i(3);
        y_start_i = line_i(4);
        x_end_i = line_i(5);
        y_end_i = line_i(6);
        
        % Find intersections with all other lines
        line_intersections = [];
        
        % 1. Check intersections with all other lines
        for j = 1:size(lines, 1)
            % Skip if it's the same line
            if acute_line_idx == j
                continue;
            end
            
            % Get parameters of the other line
            line_j = lines(j,:);
            m2 = line_j(1);
            c2 = line_j(2);
            x_start_j = line_j(3);
            y_start_j = line_j(4);
            x_end_j = line_j(5);
            y_end_j = line_j(6);
            
            % Check if lines are nearly parallel
            if abs(m1 - m2) > 1e-6
                % Calculate intersection
                x_intersect = (c2 - c1) / (m1 - m2);
                y_intersect = m1 * x_intersect + c1;
                
                % Check if intersection point is on both line segments
                on_segment_i = (x_intersect >= min(x_start_i, x_end_i) - 1e-6) && ...
                               (x_intersect <= max(x_start_i, x_end_i) + 1e-6) && ...
                               (y_intersect >= min(y_start_i, y_end_i) - 1e-6) && ...
                               (y_intersect <= max(y_start_i, y_end_i) + 1e-6);
                               
                on_segment_j = (x_intersect >= min(x_start_j, x_end_j) - 1e-6) && ...
                               (x_intersect <= max(x_start_j, x_end_j) + 1e-6) && ...
                               (y_intersect >= min(y_start_j, y_end_j) - 1e-6) && ...
                               (y_intersect <= max(y_start_j, y_end_j) + 1e-6);
                
                % If on both segments, add to line_intersections
                if on_segment_i && on_segment_j
                    line_intersections = [line_intersections; x_intersect, y_intersect, j];
                end
            end
        end
        
        % 2. Find intersections with the polygon boundary
        for j = 1:length(x_shape)-1
            % Boundary segment endpoints
            boundary_x1 = x_shape(j);
            boundary_y1 = y_shape(j);
            boundary_x2 = x_shape(j+1);
            boundary_y2 = y_shape(j+1);
            
            % Calculate intersection with boundary segment
            [xi, yi, is_intersect] = line_intersection([x_start_i, y_start_i], [x_end_i, y_end_i], ...
                                                     [boundary_x1, boundary_y1], [boundary_x2, boundary_y2]);
            
            % If valid intersection exists
            if is_intersect
                line_intersections = [line_intersections; xi, yi, -j]; % Negative j to indicate polygon boundary
            end
        end
        
        % 3. Sort intersections by distance from boundary point
        distances = zeros(size(line_intersections, 1), 1);
        for j = 1:size(line_intersections, 1)
            distances(j) = sqrt((line_intersections(j,1) - x_start_i)^2 + (line_intersections(j,2) - y_start_i)^2);
        end
        
        [~, sorted_idx] = sort(distances);
        sorted_intersections = line_intersections(sorted_idx, :);
        
        % 4. Add sorted intersections to all_points
        for j = 1:size(sorted_intersections, 1)
            x_int = sorted_intersections(j, 1);
            y_int = sorted_intersections(j, 2);
            line_idx = sorted_intersections(j, 3);
            
            % Check if point is already in all_points (avoid duplicates)
            if ~any(all(abs(all_points(:,1:2) - [x_int, y_int]) < 1e-6, 2))
                % Determine if the point is inside the bounding box
                in_bbox = (x_int >= min_x_extended - 1e-6) && (x_int <= max_x_extended + 1e-6) && ...
                          (y_int >= min_y - 1e-6) && (y_int <= max_y + 1e-6);
                
                if in_bbox
                    % Check if point is inside the polygon
                    is_inside = inpolygon(x_int, y_int, x_shape, y_shape);
                    
                    % Assign even-numbered label to intersection point
                    point_label = current_label;
                    current_label = current_label + 2; % Increment by 2 to keep even numbers
                    
                    all_points = [all_points; x_int, y_int];
                    all_intersection_info = [all_intersection_info; boundary_idx, line_idx, 1, point_label]; % 1 indicates intersection
                    
                    if is_inside
                        disp(['  Inside intersection: (', num2str(x_int), ', ', num2str(y_int), '), Label: ', num2str(point_label)]);
                        scatter(x_int, y_int, 50, 'b', 'filled');
                    else
                        disp(['  Outside intersection: (', num2str(x_int), ', ', num2str(y_int), '), Label: ', num2str(point_label)]);
                        scatter(x_int, y_int, 50, 'g', 'filled');
                    end
                    
                    % Label the intersection point
                    text(x_int, y_int+0.2, num2str(point_label), 'FontSize', 10, 'FontWeight', 'bold');
                end
            end
        end
    end
    
    % Now do the same for obtuse angle line (magenta line)
    obtuse_line_idx = find(line_from_points(:,1) == boundary_idx & line_from_points(:,2) == 1);
    
    if ~isempty(obtuse_line_idx)
        % Get line parameters
        line_i = lines(obtuse_line_idx,:);
        m1 = line_i(1);
        c1 = line_i(2);
        x_start_i = line_i(3);
        y_start_i = line_i(4);
        x_end_i = line_i(5);
        y_end_i = line_i(6);
        
        % Find intersections with all other lines
        line_intersections = [];
        
        % 1. Check intersections with all other lines
        for j = 1:size(lines, 1)
            % Skip if it's the same line
            if obtuse_line_idx == j
                continue;
            end
            
            % Get parameters of the other line
            line_j = lines(j,:);
            m2 = line_j(1);
            c2 = line_j(2);
            x_start_j = line_j(3);
            y_start_j = line_j(4);
            x_end_j = line_j(5);
            y_end_j = line_j(6);
            
            % Check if lines are nearly parallel
            if abs(m1 - m2) > 1e-6
                % Calculate intersection
                x_intersect = (c2 - c1) / (m1 - m2);
                y_intersect = m1 * x_intersect + c1;
                
                % Check if intersection point is on both line segments
                on_segment_i = (x_intersect >= min(x_start_i, x_end_i) - 1e-6) && ...
                               (x_intersect <= max(x_start_i, x_end_i) + 1e-6) && ...
                               (y_intersect >= min(y_start_i, y_end_i) - 1e-6) && ...
                               (y_intersect <= max(y_start_i, y_end_i) + 1e-6);
                               
                on_segment_j = (x_intersect >= min(x_start_j, x_end_j) - 1e-6) && ...
                               (x_intersect <= max(x_start_j, x_end_j) + 1e-6) && ...
                               (y_intersect >= min(y_start_j, y_end_j) - 1e-6) && ...
                               (y_intersect <= max(y_start_j, y_end_j) + 1e-6);
                
                % If on both segments, add to line_intersections
                if on_segment_i && on_segment_j
                    line_intersections = [line_intersections; x_intersect, y_intersect, j];
                end
            end
        end
        
        % 2. Find intersections with the polygon boundary
        for j = 1:length(x_shape)-1
            % Boundary segment endpoints
            boundary_x1 = x_shape(j);
            boundary_y1 = y_shape(j);
            boundary_x2 = x_shape(j+1);
            boundary_y2 = y_shape(j+1);
            
            % Calculate intersection with boundary segment
            [xi, yi, is_intersect] = line_intersection([x_start_i, y_start_i], [x_end_i, y_end_i], ...
                                                     [boundary_x1, boundary_y1], [boundary_x2, boundary_y2]);
            
            % If valid intersection exists
            if is_intersect
                line_intersections = [line_intersections; xi, yi, -j]; % Negative j to indicate polygon boundary
            end
        end
        
        % 3. Sort intersections by distance from boundary point
        distances = zeros(size(line_intersections, 1), 1);
        for j = 1:size(line_intersections, 1)
            distances(j) = sqrt((line_intersections(j,1) - x_start_i)^2 + (line_intersections(j,2) - y_start_i)^2);
        end
        
        [~, sorted_idx] = sort(distances);
        sorted_intersections = line_intersections(sorted_idx, :);
        
        % 4. Add sorted intersections to all_points
        for j = 1:size(sorted_intersections, 1)
            x_int = sorted_intersections(j, 1);
            y_int = sorted_intersections(j, 2);
            line_idx = sorted_intersections(j, 3);
            
            % Check if point is already in all_points (avoid duplicates)
            if ~any(all(abs(all_points(:,1:2) - [x_int, y_int]) < 1e-6, 2))
                % Determine if the point is inside the bounding box
                in_bbox = (x_int >= min_x_extended - 1e-6) && (x_int <= max_x_extended + 1e-6) && ...
                          (y_int >= min_y - 1e-6) && (y_int <= max_y + 1e-6);
                
                if in_bbox
                    % Check if point is inside the polygon
                    is_inside = inpolygon(x_int, y_int, x_shape, y_shape);
                    
                    % Assign even-numbered label to intersection point
                    point_label = current_label;
                    current_label = current_label + 2; % Increment by 2 to keep even numbers
                    
                    all_points = [all_points; x_int, y_int];
                    all_intersection_info = [all_intersection_info; boundary_idx, line_idx, 1, point_label]; % 1 indicates intersection
                    
                    if is_inside
                        disp(['  Inside intersection: (', num2str(x_int), ', ', num2str(y_int), '), Label: ', num2str(point_label)]);
                        scatter(x_int, y_int, 50, 'b', 'filled');
                    else
                        disp(['  Outside intersection: (', num2str(x_int), ', ', num2str(y_int), '), Label: ', num2str(point_label)]);
                        scatter(x_int, y_int, 50, 'g', 'filled');
                    end
                    
                    % Label the intersection point
                    text(x_int, y_int+0.2, num2str(point_label), 'FontSize', 10, 'FontWeight', 'bold');
                end
            end
        end
    end
end

% Display the total number of points collected
disp(['Total number of points collected: ', num2str(size(all_points, 1))]);

%% Step 7: Prepare for Connectivity Analysis - Use ALL points inside the bounding box
% Get labels for all points
all_labels = all_intersection_info(:, 4);

% Get all intersection points
intersection_points = all_points(all_intersection_info(:, 3) == 1, :);
intersection_labels = all_labels(all_intersection_info(:, 3) == 1);

% Get boundary points
boundary_points = all_points(all_intersection_info(:, 3) == 0, :);
boundary_labels = all_labels(all_intersection_info(:, 3) == 0);

%% Step 8: Build Connectivity Lists
% Create a mapping of all nodes to their labels
node_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
for i = 1:size(all_points, 1)
    node_map(i) = all_labels(i);
end
 
% Process Delaunay triangulation for node connectivity
dt = delaunayTriangulation(all_points(:,1), all_points(:,2));
mesh_edges = edges(dt);

% Initialize an adjacency matrix for node connectivity
num_all_nodes = max(all_labels);
adjacency_matrix = zeros(num_all_nodes, num_all_nodes);

% Add connections from mesh edges
for i = 1:size(mesh_edges, 1)
    node1 = all_labels(mesh_edges(i, 1));
    node2 = all_labels(mesh_edges(i, 2));
    adjacency_matrix(node1, node2) = 1;
    adjacency_matrix(node2, node1) = 1;
end

% Create connectivity strings in the requested format
connectivity_list = cell(size(all_points, 1), 1);
for i = 1:size(all_points, 1)
    current_node = all_labels(i);
    connections = find(adjacency_matrix(current_node, :) == 1);
    
    if ~isempty(connections)
        connectivity_list{i} = [num2str(current_node)];
        for j = 1:length(connections)
            connectivity_list{i} = [connectivity_list{i}, '-', num2str(connections(j))];
        end
    else
        connectivity_list{i} = num2str(current_node); % Node with no connections
    end
end

% Display connectivity list
fprintf('\nNode Connectivity List:\n');
for i = 1:size(all_points, 1)
    fprintf('%s\n', connectivity_list{i});
end

% Store unique node connections for reference
unique_connections = [];
for i = 1:num_all_nodes
    for j = i+1:num_all_nodes
        if adjacency_matrix(i,j) == 1
            unique_connections = [unique_connections; i, j];
        end
    end
end

connectivity_array = connectivity_list;

% Count total unique connections
fprintf('\nTotal unique connections: %d\n', size(unique_connections, 1));

%% Create the connectivity matrix as requested - MODIFIED SECTION
% Format: [node1, node2] where node1, node2 are the labels of connected nodes
connectivity = [];  % Initialize empty array instead of preallocating

% Only add connections that directly come from the Delaunay triangulation
for i = 1:size(mesh_edges, 1)
    node1 = all_labels(mesh_edges(i, 1));
    node2 = all_labels(mesh_edges(i, 2));
    
    % Ensure the smaller label is first (consistent ordering)
    if node1 > node2
        temp = node1;
        node1 = node2;
        node2 = temp;
    end
    
    % Add the connection if it's not already in the list
    if isempty(connectivity) || ~any(all(connectivity == [node1, node2], 2))
        connectivity = [connectivity; node1, node2];
    end
end

% Sort by first column, then by second column for clarity
connectivity = sortrows(connectivity);

% Display the connectivity matrix
disp('Connectivity Matrix (each row represents one edge):');
disp(connectivity);

%% Step 9: Extrude 2D Mesh into 3D
num_nodes = size(all_points, 1);
nodes_base = [all_points, zeros(num_nodes, 1)];
nodes_top  = [all_points, bridge_length * ones(num_nodes, 1)];

%% Create the coords matrix with labeled indices and their coordinates
% Format: [label, x, y, z] for both base and top nodes
coords = zeros(2*num_nodes, 4);  % Pre-allocate for efficiency

% Store base nodes first (z=0)
for i = 1:num_nodes
    coords(i, :) = [all_labels(i), nodes_base(i,1), nodes_base(i,2), nodes_base(i,3)];
end

% Store top nodes (z=bridge_length)
for i = 1:num_nodes
    coords(num_nodes+i, :) = [all_labels(i), nodes_top(i,1), nodes_top(i,2), nodes_top(i,3)];
end

disp('Created the coords matrix with format: [label, x, y, z]');
disp(['Size of coords matrix: ', num2str(size(coords,1)), ' x ', num2str(size(coords,2))]);

%% Step 10: Compute Outer Boundary
k_boundary = boundary(all_points(:,1), all_points(:,2), 0.9);

%% Step 11: Plot the Extruded 3D Mesh
figure;
hold on; grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Extruded 3D Mesh with All Intersection Points');

scatter3(nodes_base(:,1), nodes_base(:,2), nodes_base(:,3), 50, 'b', 'filled');
scatter3(nodes_top(:,1), nodes_top(:,2), nodes_top(:,3), 50, 'r', 'filled');

% Label nodes in 3D view
for i = 1:num_nodes
    text(nodes_base(i,1), nodes_base(i,2), nodes_base(i,3), num2str(all_labels(i)), 'FontSize', 10, 'FontWeight', 'bold');
    text(nodes_top(i,1), nodes_top(i,2), nodes_top(i,3), num2str(all_labels(i)), 'FontSize', 10, 'FontWeight', 'bold');
end

for i = 1:length(k_boundary)-1
    plot3(nodes_base(k_boundary(i:i+1),1), nodes_base(k_boundary(i:i+1),2), nodes_base(k_boundary(i:i+1),3), 'k-', 'LineWidth', 2);
    plot3(nodes_top(k_boundary(i:i+1),1), nodes_top(k_boundary(i:i+1),2), nodes_top(k_boundary(i:i+1),3), 'k-', 'LineWidth', 2);
end

% Plot connections based on the connectivity matrix
for i = 1:size(unique_connections, 1)
    node1_idx = find(all_labels == unique_connections(i, 1));
    node2_idx = find(all_labels == unique_connections(i, 2));
    
    if ~isempty(node1_idx) && ~isempty(node2_idx)
        plot3([nodes_base(node1_idx,1), nodes_base(node2_idx,1)], ...
              [nodes_base(node1_idx,2), nodes_base(node2_idx,2)], ...
              [nodes_base(node1_idx,3), nodes_base(node2_idx,3)], 'c-', 'LineWidth', 1.5);
        
        plot3([nodes_top(node1_idx,1), nodes_top(node2_idx,1)], ...
              [nodes_top(node1_idx,2), nodes_top(node2_idx,2)], ...
              [nodes_top(node1_idx,3), nodes_top(node2_idx,3)], 'c-', 'LineWidth', 1.5);
    end
end

for i = 1:num_nodes
    plot3([nodes_base(i,1), nodes_top(i,1)], [nodes_base(i,2), nodes_top(i,2)], [nodes_base(i,3), nodes_top(i,3)], 'm-', 'LineWidth', 1.5);
end

view(3); axis equal;
hold off;

disp(['User-defined Bridge Length: ', num2str(bridge_length), ' mm']);
disp('Upper boundary nodes labeled with odd numbers (1, 3, 5, ...)');
disp('Intersection points (both internal and external) labeled with even numbers (2, 4, 6, ...)');

% Print first few rows of coords matrix as example
disp('Sample of coords matrix (first 5 rows):');
disp('  Label    X-coord    Y-coord    Z-coord');
disp(coords(1:min(5, size(coords,1)), :));

% Save connectivity information, all points, and coords matrix to .mat file
save('mesh_connectivity.mat', 'connectivity', 'connectivity_array', 'all_labels', 'upper_labels', 'unique_connections', 'all_points', 'all_intersection_info', 'coords');
disp('Connectivity data, all points, and coords matrix saved to mesh_connectivity.mat');

%% Helper function for line intersection
function [x, y, is_intersect] = line_intersection(p1, p2, p3, p4)
    % Calculate intersection of line segment p1-p2 with line segment p3-p4
    % Returns intersection point (x,y) and whether it's a valid intersection
    
    % Line 1 parameters
    A1 = p2(2) - p1(2);
    B1 = p1(1) - p2(1);
    C1 = A1 * p1(1) + B1 * p1(2);
    
    % Line 2 parameters
    A2 = p4(2) - p3(2);
    B2 = p3(1) - p4(1);
    C2 = A2 * p3(1) + B2 * p3(2);
    
    % Calculate determinant
    det = A1 * B2 - A2 * B1;
    
    % Check if lines are parallel
    if abs(det) < 1e-6
        x = NaN;
        y = NaN;
        is_intersect = false;
        return;
    end
    
    % Calculate intersection point
    x = (B2 * C1 - B1 * C2) / det;
    y = (A1 * C2 - A2 * C1) / det;
    
    % Check if intersection is within both line segments
    if x >= min(p1(1), p2(1)) - 1e-6 && x <= max(p1(1), p2(1)) + 1e-6 && ...
       y >= min(p1(2), p2(2)) - 1e-6 && y <= max(p1(2), p2(2)) + 1e-6 && ...
       x >= min(p3(1), p4(1)) - 1e-6 && x <= max(p3(1), p4(1)) + 1e-6 && ...
       y >= min(p3(2), p4(2)) - 1e-6 && y <= max(p3(2), p4(2)) + 1e-6
        is_intersect = true;
    else
        is_intersect = false;
    end
end
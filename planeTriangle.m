function [intersectionPoints, intersectionType] = trianglePlaneIntersection3D()
    % Clear any existing figures and create a new one
    clf;
    figure('Color', 'white', 'Position', [100, 100, 1000, 600]);
    
    % User input for triangle vertices
    disp('Enter the coordinates of the triangle vertices:');
    vertices = zeros(3, 3);
    for i = 1:3
        vertices(i, :) = input(sprintf('Vertex V%d [x y z]: ', i));
    end
    
    % User input for plane
    disp('Enter the normal vector of the plane:');
    planeNormal = input('[a b c]: ');
    disp('Enter a point on the plane:');
    planePoint = input('[x y z]: ');
    
    % Calculate plane constant d
    d = -dot(planeNormal, planePoint);

    % Compute signed distances of vertices from the plane
    distances = zeros(3, 1);
    for i = 1:3
        distances(i) = dot(planeNormal(:)', vertices(i, :)) + d;
    end

    % Find intersecting points
    intersectionPoints = [];
    for i = 1:3
        j = mod(i, 3) + 1;  % Next vertex index
        
        % Line segment from vertex i to vertex j
        lineDir = vertices(j, :) - vertices(i, :);
        
        % Compute intersection point using line-plane intersection
        t = -(dot(planeNormal, vertices(i, :)) + d) / dot(planeNormal, lineDir);
        
        if t >= 0 && t <= 1
            intersectionPoint = vertices(i, :) + t * lineDir;
            intersectionPoints = [intersectionPoints; intersectionPoint];
        end
    end

    % Determine intersection type
    if isempty(intersectionPoints)
        intersectionType = 'No Intersection';
    elseif size(intersectionPoints, 1) == 1
        intersectionType = 'Intersection at Vertex';
    elseif size(intersectionPoints, 1) == 2
        intersectionType = 'Intersection along Edge';
    else
        intersectionType = 'Complex Intersection';
    end

    % 3D Visualization
    subplot(1, 2, 1);
    hold on;
    
    % Plot triangle
    trisurf([1 2 3], vertices(:,1), vertices(:,2), vertices(:,3), ...
        'FaceColor', 'blue', 'FaceAlpha', 0.3, 'EdgeColor', 'black');
    
    % Plot vertices with labels
    for i = 1:3
        plot3(vertices(i,1), vertices(i,2), vertices(i,3), 'ro', 'MarkerFaceColor', 'red', 'MarkerSize', 10);
        text(vertices(i,1), vertices(i,2), vertices(i,3), sprintf(' V%d', i), 'FontSize', 12);
    end
    
    % Plot plane
    [X, Y] = meshgrid(linspace(min(vertices(:,1))-1, max(vertices(:,1))+1, 20));
    Z = (-planeNormal(1)*X - planeNormal(2)*Y - d) / planeNormal(3);
    surf(X, Y, Z, 'FaceColor', 'red', 'FaceAlpha', 0.2);
    
    % Plot intersection points
    if ~isempty(intersectionPoints)
        plot3(intersectionPoints(:,1), intersectionPoints(:,2), intersectionPoints(:,3), ...
            'go', 'MarkerFaceColor', 'green', 'MarkerSize', 12);
        
        % Plot intersection line if exists
        if size(intersectionPoints, 1) == 2
            plot3(intersectionPoints(:,1), intersectionPoints(:,2), intersectionPoints(:,3), ...
                'g-', 'LineWidth', 3);
        end
    end
    
    title('3D Intersection Visualization');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid on;
    view(3);
    axis tight;
    
    % Signed Distances Subplot
    subplot(1, 2, 2);
    bar(distances);
    title('Signed Distances from Plane');
    xlabel('Vertex');
    ylabel('Signed Distance');
    xticklabels({'V1', 'V2', 'V3'});
    
    % Display results
    fprintf('Intersection Type: %s\n', intersectionType);
    if ~isempty(intersectionPoints)
        disp('Intersection Points:');
        disp(intersectionPoints);
    end
end

% Run the function
[intersectionPts, intersectionType] = trianglePlaneIntersection3D();
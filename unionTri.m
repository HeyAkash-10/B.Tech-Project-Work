function unionShape = unionTriangles(connectivityFile, pointsFile)
    % UNIONTRIANGLES Computes and plots the union of multiple 2D triangles
    % Inputs:
    %   connectivityFile - File containing triangle connectivity (Nx3 matrix of indices)
    %   pointsFile - File containing point coordinates (Mx2 matrix)
    
    % Read files
    connectivity = readmatrix(connectivityFile);
    points = readmatrix(pointsFile);

    % Ensure correct dimensions
    if size(points,2) ~= 2
        error('Points file must have exactly two columns (X and Y coordinates).');
    end
    if size(connectivity,2) ~= 3
        error('Connectivity file must have exactly three columns (triangle vertex indices).');
    end
    
    % Initialize an empty polyshape
    combinedShape = polyshape();

    % Loop through each triangle and add to polyshape
    for i = 1:size(connectivity,1)
        indices = connectivity(i, :);  % Get the three point indices
        x = points(indices, 1);  % Extract x-coordinates
        y = points(indices, 2);  % Extract y-coordinates
        triangle = polyshape(x, y);  % Create triangle polyshape
        combinedShape = union(combinedShape, triangle);  % Merge with previous triangles
    end

    % Plot the union of all triangles
    figure;
    plot(combinedShape, 'FaceColor', 'green', 'FaceAlpha', 0.5);
    title('Union of Triangles');
    axis equal;

    % Return the final union shape
    unionShape = combinedShape;
end

% Example usage
unionTriangles('connectivityList.txt', 'PointsList .txt');

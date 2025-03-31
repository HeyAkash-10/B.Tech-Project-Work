function unionShape = unionTriangles(x1, y1, x2, y2)
    % UNIONTRIANGLES Computes and plots the union of two 2D triangles
    % Inputs:
    %   x1, y1 - Coordinates of the first triangle (1x3 vectors)
    %   x2, y2 - Coordinates of the second triangle (1x3 vectors)

    % Create polyshape objects
    triangle1 = polyshape(x1, y1);
    triangle2 = polyshape(x2, y2);

    % Compute the union of both triangles
    unionShape = union(triangle1, triangle2);

    % Plot the original triangles
    figure;
    subplot(1,2,1);
    hold on;
    plot(triangle1, 'FaceColor', 'red', 'FaceAlpha', 0.5);
    plot(triangle2, 'FaceColor', 'blue', 'FaceAlpha', 0.5);
    title('Original Triangles');
    axis equal;
    hold off;

    % Plot the union of the triangles
    subplot(1,2,2);
    plot(unionShape, 'FaceColor', 'green', 'FaceAlpha', 0.5);
    title('Union of Triangles');
    axis equal;
end

% Define coordinates for the first triangle
x1 = [1 3 2];  
y1 = [1 1 4];  

% Define coordinates for the second triangle
x2 = [2 4 3];  
y2 = [2 2 5];  

% Call the function
unionTriangles(x1, y1, x2, y2);


function printfigure(x)

global problemID a b nballs Points LabelPoints A

% ==================================================================

if ( problemID == 2 )
    
    u = [];
    v = [];
    s = [];
    r = [];
    
    xa = [];
    xa = cell2mat(x.uv);

    u(1:nballs,1) = xa(1:2:end,1);
    v(1:nballs,1) = xa(2:2:end,1);
    s = x.s;
    r = x.r;
    
    xc = [];
    yc = [];
    for i = 1:nballs
        xc(i) = ( 1 + ( s(i) - 1 ) * b^2/a^2 ) * a * u(i);
        yc(i) = s(i) * b * v(i);
    end

    centers = [xc' yc'];

    plot_ellipse_and_circles(a, b, nballs, centers, r, [0.94,0.831,0.827]);
    return
end

% ==================================================================

if ( problemID == 3  )

    fig = figure('Visible', 'on');

    % Set the figure position and size

    set(fig, 'Position', [100, 100, 500, 500]);
    hold on;

     % Set plot limits and aspect ratio

    axis([-10 10 -10 10]);   % Define simultaneamente limites de x e y
    axis square;             % Mantém os eixos com a mesma escala

    
    % Extract and plot blue points (LabelPoints == 1)

    bluePoints = Points(:,LabelPoints == 1);
    scatter(bluePoints(1,:), bluePoints(2,:), 6, 'b','o','filled');
    
    % Extract and plot red points (LabelPoints == 0)

    redPoints = Points(:,LabelPoints == 0);
    scatter(redPoints(1,:), redPoints(2,:), 6, 'r','o','filled');
    
    % Define the matrix A and vector b

    A = x.A; 
    b = x.b;
    
    % Define the transformation matrix and ellipse center

    [Q, D] = eig(A); % Q: matrix of eigenvectors; D: diagonal matrix of eigenvalues
    theta = linspace(0, 2*pi, 100); % Parametric angle
    xy = [cos(theta); sin(theta)]; % Coordinates of a unit circle
    
    % Scale the unit circle to form the ellipse

    scale = sqrt(1 ./ diag(D)); % Scaling factors for axis lengths
    ellipse_points = Q * diag(scale) * xy;
    
    % Apply translation (shift due to vector b)

    center = -0.5 * (A \ b); % Ellipse center
    ellipse_points(1, :) = ellipse_points(1, :) + center(1);
    ellipse_points(2, :) = ellipse_points(2, :) + center(2);
    
   
    
    % Plot the ellipse

    plot(ellipse_points(1, :), ellipse_points(2, :), 'k-', 'LineWidth', 3);
    hold on;

    % Highlight the ellipse center

    scatter(center(1), center(2), 200, 'k','pentagram','filled'); 
           
end

end

% ==================================================================
% Helpers (private)
% ==================================================================

function plot_ellipse_and_circles(a, b, nballs, centers, radii, circle_colors)
    % a: semi-major axis length along x-axis
    % b: semi-minor axis length along y-axis
    % nballs: number of circles
    % centers: (nballs x 2) array with centers of the circles
    % radii: radii of the circles

    % Create figure
    figure;
    hold on;
    axis equal;
    
    % Plot the ellipse
    theta = linspace(0, 2*pi, 1000);
    x_ellipse = a * cos(theta);
    y_ellipse = b * sin(theta);
    plot(x_ellipse, y_ellipse, 'k', 'LineWidth', 2);

    % Plot the circles
    for k = 1:nballs
        x_circle = centers(k, 1) + radii * cos(theta);
        y_circle = centers(k, 2) + radii * sin(theta);
        fill(x_circle, y_circle, circle_colors, 'EdgeColor', 'k', 'LineWidth', 2);
    end
    axis off
    
    % Set axis limits
    xlim([- a - radii, a + radii]);
    ylim([- b - radii, b + radii]);

    hold off;
end
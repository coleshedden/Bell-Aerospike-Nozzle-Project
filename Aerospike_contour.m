clc; clear; close all;

% Constants
gamma = 1.4;
Me = 3.8052;                     % Exit Mach number
he_linear = 449.834;             % Linear aerospike exit height in mm
Re_annular = 150;                % Annular aerospike exit radius in mm
epsilon = 8.99334;               % Isentropic area ratio A_exit / A*

% Mach number distribution
Mx = linspace(1, Me, 300);

% Prandtl-Meyer and Mach angle functions
nu = @(M) sqrt((gamma + 1)/(gamma - 1)) .* atan( sqrt((gamma - 1)*(M.^2 - 1)/(gamma + 1)) ) ...
         - atan( sqrt(M.^2 - 1) );
mu = @(M) asin(1 ./ M);

% Exit Prandtl-Meyer angle
nu_e = nu(Me);

% Initialize arrays
hx = zeros(size(Mx));   % Linear aerospike height
xx = zeros(size(Mx));   % Linear axial location
Rx = zeros(size(Mx));   % Annular aerospike radius
Xx = zeros(size(Mx));   % Annular axial location

for i = 1:length(Mx)
    M = Mx(i);
    nu_x = nu(M);
    mu_x = mu(M);
    
    % Isentropic area ratio term
    A_ratio = (2 / (gamma + 1) * (1 + (gamma - 1) / 2 * M^2))^((gamma + 1)/(2*(gamma - 1)));
    
    % Linear aerospike height and x location
    hx(i) = he_linear * (1 - (sin(nu_e - nu_x + mu_x) / epsilon) * A_ratio);
    xx(i) = (he_linear - hx(i)) / tan(nu_e - nu_x + mu_x);
    
    % Annular aerospike radius and x location
    value_inside = 1 - (sin(nu_e - nu_x + mu_x) / epsilon) * A_ratio;
    value_inside = max(value_inside, 0);  % Prevent complex numbers
    Rx(i) = Re_annular * sqrt(value_inside);
    Xx(i) = (Re_annular - Rx(i)) / tan(nu_e - nu_x + mu_x);
end


% Plotting
figure;
plot(xx, hx, 'b-', 'LineWidth', 2); hold on;
plot(xx, -hx, 'b-', 'LineWidth', 2);              % Linear mirrored
plot(Xx, Rx, 'r--', 'LineWidth', 2);
plot(Xx, -Rx, 'r--', 'LineWidth', 2);             % Annular mirrored

% Add cowl lip lines (26.882 deg from horizontal)
theta_rad = pi/2-nu_e;

% Linear upper
x1 = xx(1);
y1 = hx(1);
y0 = y1 - x1 * tan(theta_rad);
linear_upper_lip = [0, y0, 0];  % From x = 0 to start of contour
plot([0, x1], [y0, y1], 'b-', 'LineWidth', 1.5);

% Linear lower
y1 = -hx(1);
y0 = y1 - x1 * tan(-theta_rad);
plot([0, x1], [y0, y1], 'b-', 'LineWidth', 1.5);

linear_upper_contour = [linear_upper_lip; [xx(:), hx(:), zeros(length(xx),1)]];

% Annular upper
x1 = Xx(1);
y1 = Rx(1);
y0 = y1 - x1 * tan(theta_rad);
annular_upper_lip = [0, y0, 0];
plot([0, x1], [y0, y1], 'r-', 'LineWidth', 1.5);

% Annular lower
y1 = -Rx(1);
y0 = y1 - x1 * tan(-theta_rad);
plot([0, x1], [y0, y1], 'r-', 'LineWidth', 1.5);

annular_upper_contour = [annular_upper_lip; [Xx(:), Rx(:), zeros(length(Xx),1)]];

% Final plot formatting
xlabel('x (mm)');
ylabel('Height / Radius (mm)');
title('Linear vs Annular Aerospike Contours');
legend('Linear Upper', 'Linear Lower', 'Annular Upper', 'Annular Lower', ...
       'Cowl Lip Lines');
grid on;
axis equal;

% Create 3D coordinates (z = 0) for NX ASCII format
linear_contour_ascii = [xx(:), hx(:), zeros(length(xx),1)];
annular_contour_ascii = [Xx(:), Rx(:), zeros(length(Xx),1)];

% Save as ASCII .txt files
linear_file = fopen('linear_upper_contour.txt', 'w');
fprintf(linear_file, '%.6f %.6f %.6f\n', linear_upper_contour');
fclose(linear_file);

annular_file = fopen('annular_upper_contour.txt', 'w');
fprintf(annular_file, '%.6f %.6f %.6f\n', annular_upper_contour');
fclose(annular_file);

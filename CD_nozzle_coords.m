clc; clear;

%==== NOZZLE DIMENSIONS (in meters) ====
rc = 0.110;        % Chamber radius
rt = 0.050;        % Throat radius
re = 0.150;        % Exit radius
L_conv = 0.10;    % Length of converging section
L_bell = 0.300;    % Diverging section length
theta_n = deg2rad(15); % Arc exit angle [rad]
theta_e = deg2rad(5);  % Nozzle exit angle [rad]

%==== RAO ARC SECTION ====
R_throat = 0.382 * rt;
theta_arc = linspace(0, theta_n, 50);
x_arc = R_throat * sin(theta_arc);
y_arc = rt + R_throat * (1 - cos(theta_arc));

%==== PARABOLIC DIVERGING SECTION ====
x0 = R_throat * sin(theta_n);
y0 = rt + R_throat * (1 - cos(theta_n));
A = [x0^2, x0, 1;
     L_bell^2, L_bell, 1;
     2*L_bell, 1, 0];
B = [y0; re; tan(theta_e)];
coeffs = A \ B;
a = coeffs(1); b = coeffs(2); c = coeffs(3);
x_para_early = linspace(x0, x0 + 0.01, 15);
x_para_main  = linspace(x0 + 0.01, L_bell, 100);
x_para = unique([x_para_early, x_para_main]);
y_para = a * x_para.^2 + b * x_para + c;

%==== PARABOLIC CONVERGING SECTION ====
% We'll create a parabola from x = -L_conv to 0 (throat at x = 0)
% Constraints:
%   y(-L_conv) = rc
%   y(0) = rt
%   dy/dx(0) = tan(theta_c) = slope of arc at start

theta_c = deg2rad(30);  % Arbitrary chamber angle (can be adjusted)
slope_throat = tan(theta_arc(1));  % Match arc slope at throat

A_conv = [(-L_conv)^2, -L_conv, 1;
          0, 0, 1;
          0, 1, 0];
B_conv = [rc; rt; slope_throat];
coeffs_conv = A_conv \ B_conv;
a_c = coeffs_conv(1); b_c = coeffs_conv(2); c_c = coeffs_conv(3);

x_conv = linspace(-L_conv, 0, 50);
y_conv = a_c * x_conv.^2 + b_c * x_conv + c_c;

%==== COMBINE FULL PROFILE ====
x_full = [x_conv, x_arc, x_para] * 1000;  % mm
y_full = [y_conv, y_arc, y_para] * 1000;  % mm
z_full = zeros(size(x_full));

%==== EXPORT TO TXT FILE ====
fid = fopen('nozzle_full_upper_contour_mm.txt', 'w');
for i = 1:length(x_full)
    fprintf(fid, '%.6f, %.6f, %.6f\n', x_full(i), y_full(i), z_full(i));
end
fclose(fid);
disp('Export complete: nozzle_full_upper_contour_mm.txt');

%==== PLOT PROFILE ====
figure;
plot(x_full, y_full, 'b-', 'LineWidth', 2);
xlabel('Axial Distance [mm]');
ylabel('Radius [mm]');
title('Full Rao Bell Nozzle with Smooth Converging Section');
grid on; axis equal;

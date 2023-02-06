%% Wind Tunnel Lab Sanity Check
% 6 Feb 2023

% From CRT_data_parser.m, we generate a .mat file with the following data:

%   Column  1:      Row number
%   Column  2:      Angle of attack (deg)
%   Column  3:      Elevator deflection (deg)
%   Column  4:      Rudder deflection (rad)
%   Column  5:      Air density (kg/m^3)
%   Column  6:      Air speed (m/s)
%   Column  7:      Normal force (N)
%   Column  8:      Standard deviation of normal force (N)
%   Column  9:      Transverse force (N), 
%   Column 10:      Standard deviation of transverse force (N)
%   Column 11:      Axial force (N), "X"
%   Column 12:      Standard deviation of axial force (N)
%   Column 13:      Normal Moment (N-m)
%   Column 14:      Standard deviation of normal moment (N-m)
%   Column 15:      Transverse moment (N-m)
%   Column 16:      Standard deviation of transverse moment (N-m)
%   Column 17:      Axial moment (N-m)
%   Column 18:      Standard deviation of axial moment (N-m)

% TODO: hard code actual filename
% data = load("CRT_data_.mat"); % n rows, 18 col
data = ones(2,18); % 2 rows, 18 col

% Get all elevator deflections
d = data(:,3);  % Elevator deflection (deg)
d(1) = 0;

% Only include rows where d == 0
data = data(find(d==0),:);

% Label data for clarity
a   = data(:,4);  % AoA (deg)
rho = data(:,5);  % Density (kg/m^3)
v   = data(:,6);  % Air speed (m/s)
Z   = data(:,7);  % Normal Force (N)
X   = data(:,11); % Axial Force (N)
M_n = data(:,13); % Normal moment (N*m)
M_t = data(:,15); % Traverse moment (N*m)

% Ms (sting moments) is the sum of normal and traverse moments
Ms = M_n + M_t; % (N*m)

% Given values:
c    = .2129;              % chord (m)
b_le = .9;                 % Leading edge span (m)
b_te = .985;               % Trailing edge span (m)
S    = .5*c*(b_le + b_te); % Wing area (m^2)
x_cg = .016;               % x location of center of gravity (m)
z_cg = .05;                % z location of center of gravity (m)

% We want: CL vs a, CL vs CD, and CM vs a
% First, we calculate L, D, M, and q (dynamic pressure) for each alpha.

L = X.*sind(a) - Z.*cosd(a);  % Lift
D = -X.*cosd(a) - Z.*sind(a); % Drag
M = Ms + z_cg*X + x_cg*Z;     % Moment
q = .5*rho.*(v.^2);           % Dynamic pressure 

% Next, Calculate non-dimensional coefficients
CL = L./(q*S);    % Coefficient of Lift
CD = D./(q*S);    % Coefficient of Drag
CM = M./(q*S*c);  % Moment Coefficient

% Graphs
figure()
plot(a, CL)
title("CL vs a");xlabel("a (deg)");ylabel("CL");
grid on;

figure()
plot(a, CD)
title("CD vs a");xlabel("a (deg)");ylabel("CD");
grid on;

figure()
plot(a, CM)
title("CM vs a");xlabel("a (deg)");ylabel("CM");
grid on;


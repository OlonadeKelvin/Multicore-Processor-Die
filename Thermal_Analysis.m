% Thermal Analysis of a Multi-Core Processor (Corrected)
clear; clc; close all;

%% 1. Parameters
L = 0.01;               % Chip size (m)
Nx = 50; Ny = 50;       % Grid points (x and y)
rho = 2330;             % Density (kg/m^3)
cp = 700;               % Specific heat (J/kg·K)
k = 130;                % Thermal conductivity (W/m·K)
Q_val = 5e8;            % Heat generation (W/m^3)
T_sink = 25;            % Boundary temperature (°C)
t_final = 0.1;          % Total simulation time (s)

dx = L / (Nx - 1);
dy = L / (Ny - 1);
[X, Y] = meshgrid(linspace(0, L, Nx), linspace(0, L, Ny));   % X,Y: Ny x Nx

% Initial temperature
T = ones(Ny, Nx) * T_sink;

% Core masks (logical arrays of size Ny x Nx)
core1 = (X >= 0.2*L & X <= 0.4*L & Y >= 0.2*L & Y <= 0.4*L);
core2 = (X >= 0.6*L & X <= 0.8*L & Y >= 0.6*L & Y <= 0.8*L);

%% 2. Simulation parameters (explicit Euler with safe dt)
alpha = k / (rho * cp);                 % Thermal diffusivity
dt = 1e-4;                              % Time step (satisfies stability)
steps = round(t_final / dt);            % Number of time steps

% Pre‑allocate figure
figure('Position', [100, 100, 1200, 500]);

%% 3. Time stepping loop
for s = 1:steps
    t = s * dt;
    T_old = T;

    % Set heat sources according to time
    Q = zeros(Ny, Nx);
    if t < 0.05
        Q(core1) = Q_val;
    else
        Q(core2) = Q_val;
    end

    % Update interior nodes (i: y‑index, j: x‑index)
    for i = 2:Ny-1
        for j = 2:Nx-1
            % Second derivatives
            d2Tx = (T_old(i, j+1) - 2*T_old(i, j) + T_old(i, j-1)) / dx^2;
            d2Ty = (T_old(i+1, j) - 2*T_old(i, j) + T_old(i-1, j)) / dy^2;

            % Time derivative from heat equation
            dTdt = (k * (d2Tx + d2Ty) + Q(i, j)) / (rho * cp);

            % Explicit Euler update
            T(i, j) = T_old(i, j) + dTdt * dt;
        end
    end

    % Boundary conditions remain unchanged (Dirichlet, T_sink)

    %% --- 4. Visualization every 20 steps ---
    if mod(s, 20) == 0 || s == steps
    % 2D heat map with normal y‑axis orientation
    subplot(1, 2, 1);
    imagesc([0 L]*100, [0 L]*100, T);
    colorbar;
    caxis([25 30]);
    axis xy equal tight;          % <-- normal orientation, equal scaling, tight limits
    title(sprintf('2D Heat Map (t = %.3f s)', t));
    xlabel('x (cm)'); ylabel('y (cm)');

    % 3D surface plot (no change needed)
    subplot(1, 2, 2);
    surf(X*100, Y*100, T, 'EdgeColor', 'none');
    zlabel('Temperature (°C)');
    title('3D Thermal Distribution');
    zlim([25 30]);
    view(45, 30);
    drawnow;
    end
end

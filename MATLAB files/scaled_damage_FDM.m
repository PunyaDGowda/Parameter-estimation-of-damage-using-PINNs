%% Solving Wave Equation 2D with DAMAGE MODEL - with scaling between 0 and 1
clear all;
close all;
clc;

%% Domain
% Space
original_Lx = 0.3; % x-length of the steel plate [m]
original_Ly = 0.3; % y-length of the steel plate [m]
scaled_Lx = 1.0;
scaled_Ly = 1.0;
scaling_length_factor = original_Lx*scaled_Lx;
beta = scaling_length_factor;

% wave speed (calculated from material parameters and excitation frequency)
original_c = 4773.343; % [m/s]
original_rho = 7900;
original_tsim=1.2*original_Lx/2/original_c;
scaled_tsim = 1.0;
scaling_time_factor = original_tsim*scaled_tsim;
alpha = scaling_time_factor;

scaled_c = original_c*(alpha/beta);

original_fex = 120000;
scaled_fex = original_fex * alpha;
scaled_T = 1/scaled_fex; % scaled time period

cycles=5;
scaled_dt = scaled_T*scaled_tsim/40;
original_Vamp = 10e-9; %1; %10e-9;
scaled_Vamp = original_Vamp*scaling_length_factor;

t = 0:scaled_dt:scaled_tsim;
nt = length(t);
excitation = @(s) scaled_Vamp*sin(2*pi*scaled_fex*s).*sin(pi*scaled_fex/cycles*s).^2.*(s < scaled_T*cycles);

% Space grid. dx,dy have to respect CFL, i.e. need c*dt/dx<0.5, i.e. dx > c*dt/0.5
scaled_dx = scaled_c*scaled_dt/0.5; %c*dt/0.5; % [m]
scaled_dy = scaled_dx; % [m]
x = 0:scaled_dx:scaled_Lx; % [m]
y = 0:scaled_dy:scaled_Lx; % [m]
nx = length(x);
ny = length(y);

excitation_point = scaled_Lx/2; % [m] center of the plate
excitation_index = ceil(excitation_point/scaled_dx);

%% Damage parameters
A = 0.9;
sigma_x = 0.03; % ← made these much smaller; our damages are supposed to be quite small
sigma_y = 0.03;
x_0 = 0.4; % x-coordinate of damage location [m]
y_0 = 0.3; % y-coordinate of damage location [m]

% Damage calculation
[x_d,y_d] = meshgrid(x,y);
x_part = (x_d-x_0).^2/(2*sigma_x.^2);
y_part = (y_d-y_0).^2/(2*sigma_y.^2);
damage = A*exp(-beta.^2 * (x_part + y_part));
mesh(x_d,y_d,damage)

%% Initial Conditions 
w = zeros(nx,ny,nt);
% Not defined ← sure, we have initial condition zero for the state w(:,:,2) and the derivative (i.e. w(:,:,1) = 0)
%% Saving at every 10th time step
saved_dis = zeros(nx,nx,nt);
time_data = zeros(1,nt);
saved_dis(:,:,1) = 0;
saved_dis(:,:,2) = 0;
time_data(1,1) = 0;
time_data(1,2) = scaled_dt;
snap = 3;
w = zeros(nx,ny,nt);
mesh (x,y,w(:,:,1))
xlabel('x')
ylabel('y')
%% Time stepping Loop
for n=3:nt %t=0:dt:tsim
    % Reflecting Boundary Conditions (actually not necessary…)
    w(:,[1 end],n) = 0;
    w([1 end],:,n) = 0;
    
    % forward difference of w in x-direction
    dxwp = (w(3:nx,2:ny-1,n-1) - w(2:nx-1,2:ny-1,n-1))/scaled_dx;
    % forward difference of w in y-direction
    dywp = (w(2:nx-1,3:ny,n-1) - w(2:nx-1,2:ny-1,n-1))/scaled_dy;
    % backward difference of w in x-direction
    dxwm = (w(2:nx-1,2:ny-1,n-1) - w(1:nx-2,2:ny-1,n-1))/scaled_dx;
    % backward difference of w in y-direction
    dywm = (w(2:nx-1,2:ny-1,n-1) - w(2:nx-1,1:ny-2,n-1))/scaled_dy;
    % first summand of the divergence (we average the damage, since we need it at non-integert points…)
    divx = ( (1 - 0.5*(damage(3:nx,2:ny-1)+damage(2:nx-1,2:ny-1))).*dxwp - (1 - 0.5*(damage(2:nx-1,2:ny-1)+damage(1:nx-2,2:ny-1))).*dxwm)/scaled_dx;
    % second summand of the divergence
    divy = ( (1 - 0.5*(damage(2:nx-1,3:ny)+damage(2:nx-1,2:ny-1))).*dywp - (1 - 0.5*(damage(2:nx-1,2:ny-1)+damage(2:nx-1,1:ny-2))).*dywm)/scaled_dy;
    w(2:nx-1,2:ny-1,n) = 2*w(2:nx-1,2:ny-1,n-1) - w(2:nx-1,2:ny-1,n-2) + scaled_c^2*scaled_dt^2*(divx + divy);
    
    % add the force terminal_size
    w(excitation_index,excitation_index,n) = w(excitation_index,excitation_index,n) + (scaled_dt.^2 * alpha.^2 / original_rho)*excitation(n*scaled_dt);
%         
    % Visualize at selected steps
    subplot(2,1,1)
    imagesc(x,y,w(:,:,n)');
    colorbar;
%     caxis([-1e-4 1e-4])
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('u in meters')
    title(sprintf('t = %0.10f',(n-1)*scaled_dt));
    
    subplot(2,1,2);
    mesh(x,y,w(:,:,n)) % ·
    colorbar;
%     caxis([-1e-4 1e-4])
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('u [m]')
%     axis([0 scaled_Lx 0 scaled_Ly -1e-4 1e-4]);
    pause(0)


%     if mod(n,1) == 0
%         saved_dis(:,:,snap) = w(:,:,n);
%         time_data(1,snap) = n*scaled_dt;
%         snap = snap + 1;
%    end

    saved_dis(:,:,snap) = w(:,:,n);
    time_data(1,snap) = (n-1)*scaled_dt;
    snap = snap + 1;
end

% Normalizing the displacement field with [0, 1]
min_w = min(min(min(w)));
max_w = max(max(max(w)));
scaled_w =  (w - min_w )/(max_w - min_w );

save wave2d_damage_scaled.mat scaled_w saved_dis x y time_data

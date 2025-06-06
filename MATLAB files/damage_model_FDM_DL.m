%% Solving Wave Equation 2D with DAMAGE MODEL - GOOD WORKING for REAL DATA for correct FDM
clear all;
close all;
clc;

%% Domain
% Space
Lx = 0.3; % x-length of the steel plate [m]
Ly = 0.3; % y-length of the steel plate [m]

% wave speed (calculated from material parameters and excitation frequency)
c = 4773.343; % [m/s]
rho = 7900;
%% Time and excitation
fex = 120000; % frequency in Hz used in excitation
Tex = 1/fex; % time period of the excitation [s]
cycles = 5;  % 5 full periods make up the sine burst excitation
dt =  Tex/40; % [s] (want at least 20 samples per period)
tsim = 1.2*Lx/2/c; % total simulation time [s] (slightly larger time until wave from center hits boundary, which is after Lx/2/c seconds)
t = 0:dt:tsim;
nt = length(t);
Vamp = 1e-8; % amplitude of force excitation
excitation = @(s) Vamp*sin(2*pi*fex*s).*sin(pi*fex/cycles*s).^2.*(s < Tex*cycles);

% Space grid. dx,dy have to respect CFL, i.e. need c*dt/dx<0.5, i.e. dx > c*dt/0.5
dx = c*dt/0.5;%c*dt/0.5; % [m]
dy = dx; % [m]
x = 0:dx:Lx; % [m]
y = 0:dy:Ly; % [m]
nx = length(x);
ny = length(y);

excitation_point = Lx/2; % [m] center of the plate
excitation_index = ceil(excitation_point/dx);

%% Damage parameters
A = 0.95;
sigma_x = 0.03; % ← made these much smaller; our damages are supposed to be quite small
sigma_y = 0.03;
x_0 = 0.2; % x-coordinate of damage location [m]
y_0 = 0.2; % y-coordinate of damage location [m]

% Damage calculation
[x_d,y_d] = meshgrid(x,y);
x_d_r = reshape(x_d,[151*151,1]);
y_d_r = reshape(y_d, [151*151,1]);
x_part = (x_d-x_0).^2/(2*sigma_x.^2);
y_part = (y_d-y_0).^2/(2*sigma_y.^2);
damage = A*exp(-(x_part + y_part));
mesh(damage)
%-------------------------------------
x_part_r = (x_d_r-x_0).^2/(2*sigma_x.^2);
y_part_r = (y_d_r-y_0).^2/(2*sigma_y.^2);
damage_r = A*exp(-(x_part_r + y_part_r));
%------------------------------------

% x_s = 0.15;
% y_s = 0.15;
% xx_s = x_d(2:end-1, 2:end-1);
% yy_s = y_d(2:end-1, 2:end-1);
% x_part_s = (xx_s-x_s).^2;
% y_part_s = (yy_s-y_s).^2;
% source_location = exp(-(x_part_s + y_part_s));
% %mesh(source_location)


%excitation = @(s) Vamp*sin(2*pi*fex*s).*sin(pi*fex/cycles*s).^2.*source_location *(s < Tex*cycles);
%% Saving at every 10th time step
saved_dis = zeros(nx,nx,nt);
time_data = zeros(1,nt);
saved_dis(:,:,1) = 0;
saved_dis(:,:,2) = 0;
time_data(1,1) = 0;
time_data(1,2) = dt;
snap = 3;
w = zeros(nx,ny,nt);
mesh (x,y,w(:,:,1))
xlabel('x')
ylabel('y')
%% Initial Conditions 
w = zeros(nx,ny,nt);
% Not defined ← sure, we have initial condition zero for the state w(:,:,2) and the derivative (i.e. w(:,:,1) = 0)

%% Time stepping Loop
for n=3:nt %t=0:dt:tsim
    % Reflecting Boundary Conditions (actually not necessary…)
    w(:,[1 end],n) = 0;
    w([1 end],:,n) = 0;
    
    % forward difference of w in x-direction
    dxwp = (w(3:nx,2:ny-1,n-1) - w(2:nx-1,2:ny-1,n-1))/dx;
    % forward difference of w in y-direction
    dywp = (w(2:nx-1,3:ny,n-1) - w(2:nx-1,2:ny-1,n-1))/dy;
    % backward difference of w in x-direction
    dxwm = (w(2:nx-1,2:ny-1,n-1) - w(1:nx-2,2:ny-1,n-1))/dx;
    % backward difference of w in y-direction
    dywm = (w(2:nx-1,2:ny-1,n-1) - w(2:nx-1,1:ny-2,n-1))/dy;
    % first summand of the divergence (we average the damage, since we need it at non-integert points…)
    divx = ( (1 - 0.5*(damage(3:nx,2:ny-1)+damage(2:nx-1,2:ny-1))).*dxwp - (1 - 0.5*(damage(2:nx-1,2:ny-1)+damage(1:nx-2,2:ny-1))).*dxwm)/dx;
    % second summand of the divergence
    divy = ( (1 - 0.5*(damage(2:nx-1,3:ny)+damage(2:nx-1,2:ny-1))).*dywp - (1 - 0.5*(damage(2:nx-1,2:ny-1)+damage(2:nx-1,1:ny-2))).*dywm)/dy;

    w(2:nx-1,2:ny-1,n) = 2*w(2:nx-1,2:ny-1,n-1) - w(2:nx-1,2:ny-1,n-2) + c^2*dt^2*(divx + divy);
    
    % add the force terminal_size
    w(excitation_index,excitation_index,n) = w(excitation_index,excitation_index,n) + (dt.^2/rho)*excitation(n*dt);
    %w(76,1,n) = w(76,1,n) + (dt.^2/rho)*excitation(n*dt);
%     w(2:nx-1,2:ny-1,n) = w(2:nx-1,2:ny-1,n) + (dt.^2/rho)*excitation(2:nx-1,2:ny-1,n*dt);
%         
    % Visualize at selected steps
    subplot(2,1,1)
    imagesc(x,y,w(:,:,n)');
    colorbar;
%     caxis([-1e-22 1e-22])
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('u in meters')
    title(sprintf('t = %0.10f s',n*dt));
    
    subplot(2,1,2);
    mesh(x,y,w(:,:,n)) % ·
    colorbar;
%     caxis([-1e-22 1e-22])
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('u [m]')
%     axis([0 Lx 0 Ly -1e-22 1e-22]);
    pause(0)

% 
%     if mod(n,1) == 0
%         saved_dis(:,:,snap) = w(:,:,n);
%         time_data(1,snap) = n*dt;
%         snap = snap + 1;
%     end

    saved_dis(:,:,snap) = w(:,:,n);
    time_data(1,snap) = (n-1)*dt;
    snap = snap + 1;
    
 end 
% 
 save wave2d_damage_correct_Vamp.mat saved_dis x y time_data 

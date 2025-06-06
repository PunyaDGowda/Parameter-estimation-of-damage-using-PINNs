clear all;
close all;
clc;

original_length = 0.3;  % meters
original_wave_speed = 4773.343 ; % m/s
original_fex = 120000.0 ; % Hz
n = 5;

scaling_factor_length = 1 / original_length;
scaling_factor_time = 1 / (1 / original_fex);

scaled_wave_speed = original_wave_speed * scaling_factor_time/  scaling_factor_length ;
scaled_fex = original_fex / scaling_factor_time;

T = 1 / scaled_fex;
dt = T / 40;
tsim = 1.2 * original_length / (2 * original_wave_speed);
scaled_tsim = tsim * scaling_factor_time;

Vamp = 10e-9 * scaling_factor_length;
num_steps = ceil(scaled_tsim / dt);
UU = zeros(1,500);

for i=1:200
    t=i*dt;
    UU(i)=Vamp * sin(2 * pi * scaled_fex * t) * sin(pi * scaled_fex / n * t)^2;
end

tt = 1:500;
plot(tt, UU)
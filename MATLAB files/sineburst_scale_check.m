clear all
close all
clc
%% SINEBURST EXCITATION
original_Lx = 0.3;
scaled_Lx = 1.0;
scaling_length_factor = original_Lx*scaled_Lx;
beta = scaling_length_factor;


original_c = 4773.343; % [m/s]
original_tsim=1.2*original_Lx/2/original_c;
scaled_tsim = 1.0;
scaling_time_factor = original_tsim*scaled_tsim;
alpha = scaling_time_factor;

scaled_c = original_c*(alpha/beta);

original_fex = 120000;
scaled_fex = original_fex * alpha;
scaled_T = 1/scaled_fex; % scaled time period

n=5;
scaled_dt=scaled_T*scaled_tsim/40;
original_Vamp = 10e-9;
scaled_Vamp = original_Vamp*original_Vamp;
UU = zeros(1,500);
for i=1:200
    t=i*scaled_dt;
    UU(i)=original_Vamp*sin(2*pi*scaled_fex*t)*sin(pi*scaled_fex/n*t)^2;
end

t = 1:500;
plot(t, UU)

UU_min = min(UU);
UU_max = min(UU);
scaled_UU = (UU - UU_min)/(UU_max - UU_min);
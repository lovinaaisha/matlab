%% MAIN FILE TO SOLVE AN RBC MODEL
% Lovina Putri
% Columbia University
% Macroeconomic Analysis II

clear all; close all; clc;

%% ---------- DYNARE ---------- %%
% Model parameters:
global phi alpha delta rho beta
phi = 0.7;
alpha = 0.32;
delta = 0.025;
rho = 0.95;
beta = 0.5;
save('param.mat','phi','alpha','delta','rho','beta'); % save in `param.mat'

%% IRFs for a TFP shock , as a percent deviation from steady -state
% Call Dynare to run our .mod file
dynare q3.mod

% Load the results from Dynare
load('simulation.mat')

% Plot the IRF for y and k:
irf_k = zeros(1,100); % IRFs for k_t (not k_{t+1})
irf_k(2:end)=oo_.irfs.k_eps_z(1:end-1);

% Create a tiled layout
tiledlayout(3,3);

% Output IRF
nexttile;
plot(100 * oo_.irfs.y_eps_z, 'LineWidth', 2);
title('Output (y)');

% Consumption IRF
nexttile;
plot(100 * oo_.irfs.c_eps_z, 'LineWidth', 2);
title('Consumption (c)');

% Investment IRF
nexttile;
plot(100 * oo_.irfs.I_eps_z, 'LineWidth', 2);
title('Investment (i)');

% Hours Worked IRF
nexttile;
plot(100 * oo_.irfs.l_eps_z, 'LineWidth', 2);
title('Hours Worked (l)');

% Capital IRF
nexttile;
plot(100 * irf_k, 'LineWidth', 2);
title('Capital (k)');

% Real Wage IRF
nexttile;
plot(100 * oo_.irfs.w_eps_z, 'LineWidth', 2);
title('Real Wage (w)');

% Return to Capital IRF
nexttile;
plot(100 * oo_.irfs.r_eps_z, 'LineWidth', 2);
title('Return to Capital (r)');
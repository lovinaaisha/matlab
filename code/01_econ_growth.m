%% HW 1 
% Lovina Putri
% Columbia University
% Macroeconomic Analysis II
% January 28, 2025

%% -----------------------------Problem 2--------------------------------%%
%% 2a. PDFs of log(C)
% C is lognormally distributed
% C = lognormal(mu, sigma^2)
% mu = E(log(C))
% sigma = sd(log(C))

    % log(C) has normal distribution with mean mu and variance sigma^2
    % C is in the range of (0, infty)

    % Property 1: E(C)=exp(mu+1/2*sigma^2)
        % That is: log(E(C))= mu + 1/2*sigma^2
        % Rearranging: mu = log(E(C)) - 1/2*sigma^2

% Insert the parameter
C_mean = 5; % this is E(C)
sigma_A = 0.25;
sigma_B = 1;

% NTF (mu_A, mu_B)
mu_A = log(C_mean) - 1/2*sigma_A^2;
mu_B = log(C_mean) - 1/2*sigma_B^2;

% Normal probability density for log(C):
x = -5:0.01:5;
y_A = normpdf(x, mu_A, sigma_A); 
y_B = normpdf(x, mu_B, sigma_B);

% Preferences for plots (LaTex format):
set(0, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');

% Plot
fig1 = figure(1);
plot(x, y_A, 'b', 'LineWidth', 2);
hold on
plot(x, y_B, 'r', 'linewidth', 2);
xlabel('$\log(C)$');
grid on
xticks(-5:1:5);
legend('Country A', 'Country B', 'Location','best');
set(gcf, 'color', 'w');
exportgraphics(fig1,'fig1.png');

%% 2b. PDFs of (C)

% Lognormal probability density for (C):
x = 0:0.01:20;
y_A = lognpdf(x, mu_A, sigma_A); 
y_B = lognpdf(x, mu_B, sigma_B);

% Plot
fig2 = figure(2);
plot(x, y_A, 'b', 'LineWidth', 2);
hold on
plot(x, y_B, 'r', 'linewidth', 2);
xlabel('$C$');
grid on
xticks(0:2:20);
legend('Country A', 'Country B', 'Location','best');
set(gcf, 'color', 'w');
exportgraphics(fig2,'fig2.png');

%% -----------------------------Problem 4--------------------------------%%
clear all; close all; clc
%% 4b. log(Y_t) = t*log(1+g) + log(1+eps_t)

% Parameter values:
k=0.2; f=0.2; g=0.08;

% Random number generation:
rng(450960473); 
x=rand(100,1); % 100 uniform random numbers in [0,1]
eps(x<=.5)= +k; % heads
eps(x>.5)=-f; % tails

% Plot:
t = 1:1:100;
log_Y = t*log(1+g) + log(1+eps);

fig3 = figure(3);
plot(t, log_Y, 'b', 'LineWidth', 1.75);
ylabel('log(Yt)');
xlabel('t');
grid on
set(gcf, 'color', 'w');
exportgraphics(fig3, 'fig3.png');

%% 4c. Stronger Recessions
f = f+.5; % higher value
eps(x>.5) = -f;

% Plot
t = 1:1:100;
log_Y_rec = t*log(1+g) + log(1+eps);

fig4 = figure(4);
plot(t, log_Y, 'b', 'LineWidth', 2);
hold on
plot(t, log_Y_rec, 'r', 'LineWidth', 2);
ylabel('log(Yt)');
xlabel('t');
grid on
legend('Normal Recessions', 'Stronger Recessions', 'Location','best');
exportgraphics(fig4, 'fig4.png');

%% 4d. Faster Growth Rate
g = 0.15;
f = 0.2; eps(x>.5)=-f;

% Plot:
t = 1:1:100;
log_Y_g = t*log(1+g) + log(1+eps);

fig5 = figure(5);
plot(t, log_Y, 'b', 'LineWidth', 2);
hold on
plot(t, log_Y_g, 'r', 'LineWidth', 2);
ylabel('log(Yt)');
xlabel('t');
grid on
legend('Normal Growth Rate', 'Faster Growth Rate', 'Location', 'best');
set(gcf, 'color', 'w');
exportgraphics(fig5, 'fig5.png');

%% 4e. Persistent Output
% Generate persistent cycles:
e = NaN(1, 100);
e(1)=eps(1);
for i = 2:100
    if e(i-1)==+k
        e(x<=.9)=+k;
        e(x>.9)=-f;
    else
        e(x<=.1)=+k;
        e(x>.1)=-f;
    end
end

% Plot:
g = 0.08;
t = 1:1:100;
log_Y_pers = t*log(1+g) + log(1+e);

fig6 = figure(6);
plot(t, log_Y,'b', 'LineWidth', 2);
hold on
plot(t, log_Y_pers, 'r', 'LineWidth', 2);
ylabel('log(Yt)');
xlabel('t');
grid on
legend('Normal Output', 'Persistent Output', 'Location', 'best');
set(gcf, 'color', 'w');
exportgraphics(fig6, 'fig6.png');

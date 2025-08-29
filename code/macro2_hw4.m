%% Problem Set 4
% Lovina Putri
% Columbia University
% Macroeconomic Analysis II
% 20 February, 2025
clear all; close all; clc;

%% ---------- DATA ---------- %%
% Import the .csv FRED data:
gdp = readtable("GDPC1.csv"); % Industrial Production [seasonally adjusted]
gcec = readtable("GCEC1.csv"); % Nonfarm employees (000s) [seasonally adjusted]

% Rename the variables:
gdp.Properties.VariableNames= {'date','gdp'};
gcec.Properties.VariableNames= {'date','gcec'};

% Read the dates in date format:
gdp.date = datetime(gdp.date);
gcec.date = datetime(gcec.date);

% Merge the two datasets by date, in a matrix called 'data':
data = outerjoin(gdp,gcec,'LeftKeys','date','Rightkeys','date','mergekeys',true);

% Keep observations from 1969:1 to 2024:10
sdate=datetime('01-Jan-1969');
edate=datetime('01-Oct-2024');  
data_start=find(datenum(data.date)==datenum(sdate),1,'first');
data_end=find(datenum(data.date)==datenum(edate),1,'first');
data = data(data_start:data_end,:);

% Remove missing values:
data = rmmissing(data); % remove rows with NaN
data.Properties.VariableNames = {'date', 'gdp', 'gcec'}; % (new) variable names

% Substract each individual variable and generate a time-series vector, Y:
D=data; date = D.date; gdp = D.gdp; gcec = D.gcec;
Y = [gdp, gcec];

% Take logs of the time series:
log_gdp = log(Y(:,1));
log_gcec = log(Y(:,2));

%% ---------- PLOTS: RAW TIME SERIES ---------- %%
% Preferences for plots (LaTeX format)
set(0,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

% Log GDP:
plot(date,log_gdp,'b','LineWidth',2);
ylabel('ln(Yt)');
xlabel('Period $t$');
grid on
set(gcf,'color','w');

% Log GCEC:
close all
plot(date,log_gcec,'b','LineWidth',2);
ylabel('ln(Gt)');
xlabel('Period $t$');
grid on
set(gcf,'color','w');

%% ---------- TREND VS. CYCLE ---------- %%
% function: hpfilter(Time-series,Lambda)
%
% separates one or more time series into trend and cyclical components with the
% Hodrick-Prescott filter.
%
%	[Trend, Cycle] = hpfilter(Series, Smoothing);
%
% Outputs:
%	Trend - Trend component of Series.
%	Cycle - Cyclical component of Series.
%
% Notes:
%	The Hodrick-Prescott filter separates a time series into a "trend"
%	component and a "cyclical" component such that
%		Series = Trend + Cyclical
%
%	The reference suggests values for the smoothing parameter that depend upon
%	the periodicity of the data with:
%		Periodicity		Smoothing
%		Yearly			100
%		Quarterly		1600
%		Monthly			14400
%

% Smoothing parameter:
lambda = 1600; % quarterly

% HP filter:
[trend_gdp,cycle_gdp] = hpfilter(log_gdp,lambda);
[trend_gcec,cycle_gcec] = hpfilter(log_gcec,lambda);

% time-series vs. trend:
% Log GDP:
close all
plot(date,log_gdp,'b','LineWidth',2);
hold on
plot(date,trend_gdp,'r','LineWidth',1.5);
ylabel('Log Values');
xlabel('$t$');
grid on
set(gcf,'color','w');
legend('ln(Yt)','Trend','location','best');

% Log GCEC:
close all
plot(date,log_gcec,'b','LineWidth',2);
hold on
plot(date,trend_gcec,'r','LineWidth',1.5);
ylabel('Log Values');
xlabel('$t$');
grid on
set(gcf,'color','w');
legend('ln(Gt)','Trend','location','best');

%% ---------- DETRENDED DATA ---------- %%
gdp_detrended = log_gdp - trend_gdp;
gcec_detrended = log_gcec - trend_gcec;

% Plot:
% Log GDP
close all
plot(date,gdp_detrended,'b','LineWidth',2);
ylabel('Detrended ln(Yt)');
xlabel('$t$');
grid on
set(gcf,'color','w');

% Log GCEC
close all
plot(date,gcec_detrended,'b','LineWidth',2);
ylabel('Detrended ln(Gt)');
xlabel('$t$');
grid on
set(gcf,'color','w');

%% Find the standard deviation of a de-trended time series
% Log GDP:
gdp_detrended_std = std(gdp_detrended)

% Log GCEC:
gcec_detrended_std = std(gcec_detrended)

% Ratio:
ratio_gdptogc = gdp_detrended_std/gcec_detrended_std
ratio_gctogdp = gcec_detrended_std/gdp_detrended_std

%% Find the correlation of lagged de-trended time series
% G is a T-by-1 vector
T = length(gcec_detrended);

% Preallocate with NaNs to keep the same length
gcec_lag = nan(T,1);

% Shift each observation down by 1
gcec_lag(2:end) = gcec_detrended(1:end-1);

% G is a T-by-1 vector
T1 = length(gcec_detrended);

% Preallocate G_lead (use NaNs to keep dimensions consistent)
gcec_lead = nan(T1, 1);

% Shift the vector upward:
gcec_lead(1:end-1) = gcec_detrended(2:end);


%% Find the correlation of 2 de-trended time series
corr_detrended = corr(gdp_detrended,gcec_detrended)

% Correlation
corr_detrended_1 = corr(gdp_detrended(2:end),gcec_lag(2:end))

% Correlation
corr_detrended_2 = corr(gdp_detrended(1:end-1), gcec_lead(1:end-1))

% Plot:
close all
plot(date,gdp_detrended,'b','LineWidth',2);
hold on
plot(date,gcec_detrended,'r','LineWidth',2);
ylabel('Log Values');
xlabel('$t$');
grid on
set(gcf,'color','w');
legend('Detrended ln(Yt)','Detrended ln(Gt)','location','southwest');
% --------------------------------------------------------------------- END
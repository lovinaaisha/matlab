%% HW7
% Lovina Putri
% Columbia University
% Macroeconomic Analysis II

%% ---------- Question 2 ---------- %%
clear all; close all; clc;

% Read the data into a table (assumes the first row contains headers)
data = readtable("fred.xlsx");

% Convert 'date' column from text to datetime
data.date = datetime(data.date, 'InputFormat','dd-MMM-yyyy');

% Extract the relevant columns from the table
inflation_rate = data.inflation_rate;
tbill_rate = data.tbill_rate;

% Calculate the correlation coefficient
r = corr(tbill_rate, inflation_rate);

% Display the result
fprintf('The correlation coefficient between inflation rate and T-bill rate is: %f\n', r);

% Plotting the two series on the same line graph
figure;
plot(data.date, inflation_rate, '-b', 'LineWidth', 2);  % Plot inflation rate in blue
hold on;
plot(data.date, tbill_rate, '-r', 'LineWidth', 2);        % Plot T-bill rate in red
xlabel('Date');
ylabel('Rate');
title('Inflation Rate (%YoY) vs. 3-months T-bills yields');
legend('Inflation Rate', 'T-bills Yield', 'Location', 'best');
grid on;
hold off;


%% ---------- Question 3 (a) ---------- %%
clear
CPI3data = readtable('[path]/Macro/HW7/CPI3.csv')
M2data = readtable('[path]/Macro/HW7/M2.csv')
CPI3 = CPI3data.CPIAUCSL
M2 = M2data.M2SL
date = CPI3data.observation_date
inflation3 = (CPI3(2:end) ./ CPI3(1:end-1) - 1) * 100
M2rate = (M2(2:end) ./ M2(1:end-1) - 1) * 100
CPI3data.inflation3 = [NaN; inflation3]
M2data.M2rate = [NaN; M2rate]

%% ---------- Question 3 (b) ---------- %%
plot(date(2:end),M2rate,'b','LineWidth',2);
hold on
plot(date(2:end),inflation3,'c','LineWidth',2);
hold on
ylabel('Inflation and Money growth rate');
xlabel('Year');
grid on
set(gcf,'color','w');
legend('Money growth rate','Inflation','location','best');

corr3 = corr (M2rate, inflation3)
%% ---------- Question 3 (c) function ---------- %%
function [Trend, Cyclical] = hpfilter(Series, Smoothing)

% Step 1 - check arguments

if nargin < 1 || isempty(Series)
	error('GARCH:hpfilter:MissingInputData', ...
		'Required times series data Series missing or empty.');
end

if ~isscalar(Series) && isvector(Series) && isa(Series,'double')
	Series = Series(:);
	[NumSamples, NumSeries] = size(Series);
elseif ndims(Series) == 2 && min(size(Series)) > 1 && isa(Series,'double')
	[NumSamples, NumSeries] = size(Series);
else
	error('GARCH:hpfilter:InvalidInputArg', ...
 		'Invalid format for time series data Series. Must be a vector or matrix.');
end

if any(any(~isfinite(Series)))
	error('GARCH:hpfilter:InvalidInputData', ...
		'Cannot have infinite or missing values in data.');
end

if NumSamples < 3		% treat samples with < 3 observations as trend data only
	warning('GARCH:hpfilter:InsufficientData', ...
		'Need at least three observations. Will just return input as Trend.');
	Trend = Series;
	Cyclical = zeros(NumSamples, NumSeries);
	return
end

if nargin < 2 || isempty(Smoothing)
	warning('GARCH:hpfilter:DefaultQuarterlySmoothing', ...
		'Missing or empty Smoothing parameter set to 1600 (quarterly data).');
	Smoothing = 1600;
end

if ~isvector(Smoothing) || ~isa(Smoothing,'double')
	error('GARCH:hpfilter:InvalidInputArg', ...
 		'Invalid format for Smoothing parameter. Must be a scalar or vector.');
else
	if ~any(numel(Smoothing) == [1, NumSeries])
		error('GARCH:hpfilter:InconsistentSmoothingDimensions', ...
			'Smoothing parameter is neither a scalar nor a conformable vector.');
	end
end

if any(isnan(Smoothing))
	error('GARCH:hpfilter:InvalidSmoothing', ...
		'Must have finite or infinite smoothing parameter.');
end

if any(Smoothing < 0)
	warning('GARCH:hpfilter:NegativeSmoothing', ...
		'Negative value for smoothing parameter. Will use absolute value.');
	Smoothing = abs(Smoothing);
end

% Step 2 - run the filter with either scalar or vector Smoothing

if (numel(Smoothing) == 1) || (max(Smoothing) == min(Smoothing))	% scalar
	if numel(Smoothing) > 1
		Smoothing = Smoothing(1);
	end
	if isinf(Smoothing)			% do OLS detrending
		Trend = Series - detrend(Series);
	else
		if NumSamples == 3			% special case with 3 samples
			A = eye(NumSamples, NumSamples) + ...
				Smoothing*[ 1 -2 1; -2 4 -2; 1 -2 1 ];
		else						% general case with > 3 samples
			e = repmat([Smoothing, -4*Smoothing, (1 + 6*Smoothing), ...
				-4*Smoothing, Smoothing], NumSamples, 1);
			A = spdiags(e, -2:2, NumSamples, NumSamples);
			A(1,1) = 1 + Smoothing;
			A(1,2) = -2*Smoothing;
			A(2,1) = -2*Smoothing;
			A(2,2) = 1 + 5*Smoothing;
			A(NumSamples - 1, NumSamples - 1) = 1 + 5*Smoothing;
			A(NumSamples - 1, NumSamples) = -2*Smoothing;
			A(NumSamples, NumSamples - 1) = -2*Smoothing;
			A(NumSamples, NumSamples) = 1 + Smoothing;
		end
		Trend = A \ Series;
	end
else																% vector
	Trend = zeros(NumSamples,NumSeries);
	if NumSamples == 3						% special case with 3 samples
		for i = 1:NumSeries
			if isinf(Smoothing(i))			% do OLS detrending
				Trend(:, i) = Series(:, i) - detrend(Series(:, i));
			else
				A = eye(NumSamples, NumSamples) + ...
					Smoothing(i)*[ 1 -2 1; -2 4 -2; 1 -2 1 ];
				Trend(:, i) = A \ Series(:, i);
			end
		end
	else									% general case with > 3 samples
		for i = 1:NumSeries
			if isinf(Smoothing(i))			% do OLS detrending
				Trend(:, i) = Series(:, i) - detrend(Series(:, i));
			else
				e = repmat([Smoothing(i), -4*Smoothing(i), (1 + 6*Smoothing(i)), ...
					-4*Smoothing(i), Smoothing(i)], NumSamples, 1);
				A = spdiags(e, -2:2, NumSamples, NumSamples);
				A(1,1) = 1 + Smoothing(i);
				A(1,2) = -2*Smoothing(i);
				A(2,1) = -2*Smoothing(i);
				A(2,2) = 1 + 5*Smoothing(i);
				A(NumSamples - 1, NumSamples - 1) = 1 + 5*Smoothing(i);
				A(NumSamples - 1, NumSamples) = -2*Smoothing(i);
				A(NumSamples, NumSamples - 1) = -2*Smoothing(i);
				A(NumSamples, NumSamples) = 1 + Smoothing(i);
				Trend(:, i) = A \ Series(:, i);
			end
		end
	end
end

% Step 3 - plot results if no output arguments

if nargout == 0
	figure(gcf);
	plot(Series,'b');
	hold on
	plot(Trend,'r');
	title('\bfHodrick-Prescott Filter');
	if NumSeries == 1
		legend('Raw Data','Smoothed Trend');
	end
	hold off;
elseif nargout > 1
	Cyclical = Series - Trend;
end
end
%% ---------- Question 3 (c) ---------- %%
lambda = 100 %Annual
[trend_M2rate,cycle_M2rate] = hpfilter(M2rate,lambda)
[trend_inflation3,cycle_inflation3] = hpfilter(inflation3,lambda)
corr_trend = corr(trend_inflation3,trend_M2rate)
plot(date(2:end),M2rate,'b','LineWidth',2);
hold on
plot(date(2:end),trend_M2rate,'r','LineWidth',1);
hold on
plot(date(2:end),inflation3,'c','LineWidth',2);
hold on
plot(date(2:end),trend_inflation3,'m','LineWidth',1);
ylabel('Inflation and money growth rate');
xlabel('Year');
grid on
set(gcf,'color','w');
legend('Money growth rate','Money growth rate trend','Inflation','Inflation trend','location','best');

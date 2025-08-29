%% PROGRAM TO SOLVE AN RBC MODEL WITH DYNARE

%% Labeling block
% Endogenous variables (8) - in logs :
var
y % output
I % investment
k % capital
l % labor
z % productivity
c % consumption
w % real wage
r % real interest rate
;

% Exogenous variables (1):
varexo
eps_z % productivity shock
;

%% Parameters block
% Define parameters :
parameters phi alpha delta rho beta ; % model parameters
parameters r_ss l_ss k_ss y_ss w_ss c_ss i_ss z_ss ; % steady - state values

% Set model parameter values loading the `param .mat ' file :
load param.mat ;
set_param_value ('phi',phi); set_param_value ('alpha',alpha);
set_param_value ('delta',delta); set_param_value ('rho',rho);
set_param_value ('beta',beta);

% Set steady - state values using steady - state closed - form equations :
r_ss = 1/ beta -(1- delta );
l_ss = 1/(1+ phi /(1 - phi )/(1 - alpha )*(1 - alpha * delta / r_ss ));
k_ss = l_ss *( alpha / r_ss )^(1/(1 - alpha ));
y_ss = k_ss^alpha * l_ss^(1 - alpha );
w_ss = (1- alpha )* y_ss / l_ss ;
c_ss = y_ss - delta * k_ss ;
i_ss = delta * k_ss ;
z_ss = 1;

%% Model block
% Equilibrium conditions in logs , for log - linearization (8) :
model ;

% Euler equation :
(1- phi )/exp (c)= beta *(1 - phi )/exp (c (+1) )*(1 - delta +exp(r (+1) ));

% Labor-supply :
(phi)/ (1 - exp(l)) = exp(w)*((1-phi)/exp(c));

% Production function :
exp (y)=exp(z)*exp(k( -1))^alpha*exp(l)^(1 - alpha );

% Goods market clearing:
exp (y)= exp(c) + exp(I); 

% FOC for capital:
exp(r) = alpha*(exp(y))/(exp(k));

% FOC for labor:
exp(w) = (1-alpha)*(exp(y))/(exp(l));

% Law of motion for capital:
exp(k) = (1-delta)*exp(k(-1)) + exp(I);

% Productivity equation ( already log - linearized )
z = rho*z( -1) + eps_z ;

end ;

%% Initialization block
% Initial guess for the endogenous variables in logs (8)
initval ;
y = log( y_ss ); % steady - state log - value of output
I = log( i_ss ); % steady - state log - value of investment
k = log( k_ss ); % steady - state log - value of capital
l = log( l_ss ); % steady - state log - value of labor
z = log( z_ss ); % steady - state log - value of productivity
c = log( c_ss ); % steady - state log - value of consumption
w = log( w_ss ); % steady - state log - value of real wage
r = log( r_ss ); % steady - state log - value of real interest rate
end ;

% Compute the steady state :
steady(maxit=10000) ; % maximum number of iterations
check; % check that the steady state is indeed found
resid; % inspect the model equations residuals ( all should be 0)

% Variance of exogenous shocks (1):
shocks;
var eps_z; stderr log (1.01); 
end;

%% Solution block
% Perform the simulation of the model . Get results for specified variables
stoch_simul(order=1,irf=100,noprint ,nodisplay) y c I l k w r z;

% Store the results :
savefile = 'simulation.mat ';
save (savefile, 'oo_');


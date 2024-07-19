%% This code was written and run in Octave. No guaruntees that it will be compatible with MATLAB
%% Author Ryan Pewowaruk - Ryan Pewowaruk Research Consulting

%% Local Carotid Artery Analysis

%% In the actual analysis data is imported from CSV file.

%% Assign NTG Data for Example
Dd_NTG = 7.35;
Ds_NTG = 7.71;
Pd_NTG = 70;
Ps_NTG = 114;

%% Assign Baseline Data
Dd = 6.73;
Ds = 7.29;
Pd = 77;
Ps = 135;

kref = 0.1; % reference k of 0.1
Pref = 80; % reference pressure of 80mmHg

beta_ECM = log(Ps_NTG./Pd_NTG)/(Ds_NTG./Dd_NTG-1)-log(Pd_NTG./Pref); % Calculate beta ECM analytically
D0 = Dd_NTG./(log(Pd_NTG/Pref)./beta_ECM+1); % Calculate reference diameter analytically

%% Make array to fill in for each participant in for loop
output_carotid = zeros(length(beta_ECM), 5);

for i=1:length(beta_ECM)

% Create optimization function to solve for beta VSMC

  fx = @(param) Pref * [exp(beta_ECM(i)*(Dd(i)/D0(i) -1)) + param(1)/kref * exp(param(2)*(Dd(i)/((1-param(1))*D0(i)) -1)) ...
  exp(beta_ECM(i)*(Ds(i)/D0(i) -1)) + param(1)/kref * exp(param(2)*(Ds(i)/((1-param(1))*D0(i)) -1)) ] ...
  - [Pd(i) Ps(i)];

  % Solve optimization function
  [param, fvec, info] = fsolve(fx, [0.05 10]);

 % Assign variables into output matrix (output_cfPWV)
  % Column 1 - Dref (Reference diameter)
  % Column 2 - betaECM (ECM stiffness)
  % Column 3 - betaVSMC (VSMC stiffness)
  % Column 4 - k (VSMC Tone)
  % Column 5 - Output if model solved properly. See help fsolve for more information
output_carotid(i,:) = [D0(i), beta_ECM(i) param(2) param(1) info];

  end

%% cfPWV Analysis

%% In the actual analysis data is imported from CSV file.

%% Assign example cfPWV data (Pressure is diastolic BP)
P_Base = 77;
P_NTG = 72;
PWV_Base = 8.5;
PWV_NTG = 9.6;

%% Make array to fill in for each participant in for loop
output_cfPWV = zeros(length(P_Base), 4);

k=0.05; % Assume k=0.05 for all participants
kref = 0.1; % Assign reference k value of 0.1

for i=1:length(P_Base)

beta_ECM=2*1050*PWV_NTG(i)^2/(133.32*P_NTG(i))-log(P_NTG(i)/80); % Analytic calculation for beta ECM

% Create optimization function to solve for beta VSMC
fx = @(param) [80*(exp(beta_ECM*(param(1)-1))+k/kref*exp(param(2)*(param(1)/(1-k)-1))) ...
sqrt(133.32*80/(2*1050)*(beta_ECM*exp(beta_ECM*(param(1)-1))+k/kref*param(2)/(1-k)*exp(param(2)*(param(1)/(1-k)-1)))) ] ...
-[P_Base(i) PWV_Base(i)];

% Solve optimization function
[param, fvec, info] = fsolve(fx, [0.90 10]);

% Assign variables into output matrix (output_cfPWV)
  % Column 1 - betaECM (ECM stiffness)
  % Column 2 - betaVSMC (VSMC stiffness)
  % Column 3 - DDr (Ratio of diameter to reference diameter)
  % Column 4 - Output if model solved properly. See help fsolve for more information

output_cfPWV(i,:) = [beta_ECM param(2) param(1) info];
end






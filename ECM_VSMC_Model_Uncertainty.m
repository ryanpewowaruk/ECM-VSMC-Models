%% This code was written and run in Octave. No guaruntees that it will be compatible with MATLAB
%% Author Ryan Pewowaruk - Ryan Pewowaruk Research Consulting

% Senstivity analysis and uncertainty quantification for parameter estimation of ECM and smooth muscle
% stiffness parameters using arterial mechanics model. Both uncertainty from
% measurement error as well as modeling assumptions are assessed

%%%% Local Carotid Artery Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Measurement Uncertainty (1 Standard Deviation)
U_Dd = 0.05;
U_Ds = 0.05;
U_Pd = 5;
U_Ps = 7;
U_PWV = 0.8;

%% Assign NTG Data
Dd_NTG = repmat(7.35, [18, 1]);
Ds_NTG = repmat(7.71, [18, 1]);
Pd_NTG = repmat(70, [18, 1]);
Ps_NTG = repmat(114, [18, 1]);

%% Assign Baseline Data
Dd = repmat(6.73, [18, 1]);
Ds = repmat(7.29, [18, 1]);
Pd = repmat(77, [18, 1]);
Ps = repmat(135, [18, 1]);

%% Add in Measurement Uncertainty - One at a time method
Dd_NTG(1:2) = Dd_NTG(1:2) + [U_Dd; -U_Dd];

Dd(3:4) = Dd(3:4) + [U_Dd; -U_Dd];

Ds_NTG(5:6) = Ds_NTG(5:6) + [U_Ds; -U_Ds];
Ds(7:8) = Ds(7:8) + [U_Ds; -U_Ds];

Pd_NTG(9:10) = Pd_NTG(9:10) + [U_Pd; -U_Pd];
Pd(11:12) = Pd(11:12) + [U_Pd; -U_Pd];

Ps_NTG(13:14) = Ps_NTG(13:14) + [U_Ps; -U_Ps];
Ps(15:16) = Ps(15:16) + [U_Ps; -U_Ps];

%% Assign Reference Values
kref = 0.1; % reference k of 0.1
Pref = 80; % reference pressure of 80mmHg

beta_ECM = log(Ps_NTG./Pd_NTG)./(Ds_NTG./Dd_NTG-1)-log(Pd_NTG./Pref); % Calculate beta ECM analytically
D0 = Dd_NTG./(log(Pd_NTG/Pref)./beta_ECM+1); % Calculate reference diameter analytically

%% Make output array to fill in for each participant in for loop
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

%% Assess assumption that NTG eliminates 100% of VSMC tone (change from 100% to 80%)
 for i=18

   fx = @(param) Pref * [exp(param(3)*(Dd_NTG(i)/param(4) -1)) + 0.2*param(1)/kref * exp(param(2)*(Dd_NTG(i)/((1-0.2*param(1))*param(4)) -1)) ...
  exp(param(3)*(Ds_NTG(i)/param(4) -1)) + 0.2*param(1)/kref * exp(param(2)*(Ds_NTG(i)/((1-0.2*param(1))*param(4)) -1))  ...
  exp(param(3)*(Dd(i)/param(4) -1)) + param(1)/kref * exp(param(2)*(Dd(i)/((1-param(1))*param(4)) -1)) ...
  exp(param(3)*(Ds(i)/param(4) -1)) + param(1)/kref * exp(param(2)*(Ds(i)/((1-param(1))*param(4)) -1)) ]...
  - [Pd_NTG(i) Ps_NTG(i) Pd(i) Ps(i)];

  % Solve optimization function
  [param, fvec, info] = fsolve(fx, [0.05 10 10 6]);

 % Assign variables into output matrix (output_cfPWV)
  % Column 1 - Dref (Reference diameter)
  % Column 2 - betaECM (ECM stiffness)
  % Column 3 - betaVSMC (VSMC stiffness)
  % Column 4 - k (VSMC Tone)
  % Column 5 - Output if model solved properly. See help fsolve for more information
output_carotid(i,:) = [param(4), param(3) param(2) param(1) info];

   end

%% Calculate uncertainty for each parameter permutation
U_carotid = output_carotid(2:2:end, 1:4) - output_carotid(1:2:end, 1:4);

%% Calculate cumulative uncertainties
Uc_D_meas = sqrt(sum(U_carotid(1:4, :).^2) ); % Diameter Measurement Uncertainty
Uc_P_meas = sqrt(sum(U_carotid(5:8, :).^2) ); % Pressure Measurement Uncertainty
Uc_meas = sqrt(sum(U_carotid(1:8, :).^2) ); % Total Measurement Uncertainty
Uc_assump = U_carotid(9,:); % Model Assumptions Uncertainty

Uc_total = sqrt( Uc_meas.^2 + Uc_assump.^2); % Total Uncertainty


%%%%% cfPWV Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Make vector of cfPWV and BP data (Pressure is diastolic BP)

P_Base = repmat(76, [12,1]);
P_NTG = repmat(71.5, [12,1]);
PWV_Base = repmat(7.05, [12,1]);
PWV_NTG = repmat(7.5, [12,1]);

%% Add in Measurement Uncertainty - One at a time method
P_Base(1:2) = P_Base(1:2) + [U_Pd; -U_Pd];
P_NTG(3:4) = P_NTG(3:4) + [U_Pd; -U_Pd];

PWV_Base(5:6) = PWV_Base(5:6) + [U_PWV; -U_PWV];
PWV_NTG(7:8) = PWV_NTG(7:8) + [U_PWV; -U_PWV];

%% Make output array to fill in for each participant in for loop
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


%% Assess assumption that k=0.05 for all participatns (change from 0.05 to 0.02)
for i=10

k=0.02;

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

%% Assess assumption that NTG eliminates 100% of VSMC tone (change from 100% to 80%)
for i=12

k=0.05;
k_ntg = 0.01; %(0.01 is 80% reduction from 0.05)

beta_ECM=2*1050*PWV_NTG(i)^2/(133.32*P_NTG(i))-log(P_NTG(i)/80); % Analytic calculation for beta ECM

% Create optimization function to solve for beta VSMC
fx = @(param) [80*(exp(param(3)*(param(1)-1))+k/kref*exp(param(2)*(param(1)/(1-k)-1))) ...
sqrt(133.32*80/(2*1050)*(param(3)*exp(param(3)*(param(1)-1))+k/kref*param(2)/(1-k)*exp(param(2)*(param(1)/(1-k)-1)))) ...

80*(exp(param(3)*(param(4)-1))+k_ntg/kref*exp(param(2)*(param(4)/(1-k_ntg)-1))) ...
sqrt(133.32*80/(2*1050)*(param(3)*exp(param(3)*(param(4)-1))+k_ntg/kref*param(2)/(1-k)*exp(param(2)*(param(4)/(1-k_ntg)-1)))) ] ...

-[P_Base(i) PWV_Base(i) P_NTG(i) PWV_Base(i)];

% Solve optimization function
[param, fvec, info] = fsolve(fx, [0.90 10 10, 0.95]);

% Assign variables into output matrix (output_cfPWV)
  % Column 1 - betaECM (ECM stiffness)
  % Column 2 - betaVSMC (VSMC stiffness)
  % Column 3 - DDr (Ratio of diameter to reference diameter)
  % Column 4 - Output if model solved properly. See help fsolve for more information

output_cfPWV(i,:) = [param(3) param(2) param(1) info];
end

%% Calculate uncertainty for each parameter permutation
U_cfPWV= output_cfPWV(2:2:end, 1:2) - output_cfPWV(1:2:end, 1:2);

%% Calculate cumulative uncertainties
Ucf_P_meas = sqrt(sum(U_cfPWV(1:2, :).^2) ); % Pressure Measurement Uncertainty
Ucf_PWV_meas = sqrt(sum(U_cfPWV(3:4, :).^2) ); % PWV Measurement Uncertainty
Ucf_meas = sqrt(sum(U_cfPWV(1:4, :).^2) ); % Total Measurement Uncertainty
Ucf_assump = sqrt(sum(U_cfPWV(5:6, :).^2) ); % Model Assumption Uncertainty

Ucf_total = sqrt( Ucf_meas.^2 + Ucf_assump.^2); % Total Uncertainty

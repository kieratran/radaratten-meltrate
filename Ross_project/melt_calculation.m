% ==== Author: Kiera Tran ====
% This script generate all basal melt rate vs attenuation rate profiles simulated from D. Yang's 1.5D thermal model
% % !!! WARNING !!! : Each tmd_Na_SD and tmd_Na_TD is a vector, which can add up very fast and take up plenty of storage memory. 
% Hence, MATLAB could crash

%% Generate basal melt vs attenuation relationship:

% ---- 1. Loading inputs for the model ----
% % !!!! Matlab will crash if the dataset is too big !!! % %
load('Model_inputs.mat');
thickness = H(1:100); velocity = v(1:100); basal_slope = bs(1:100); surface_slope = ss(1:100); 
smb = smb(1:100); dTdl =dTdl(1:100); Ts =Ts(1:100);
% thickness = total ice thickness (m)
% velocity = ice velocity (m/yr) from MEaSUREs version 2
% dTdl = horizontal temperature gradient (K/yr) from E. Dawson's ISSM thermal model
% smb = surface mass balance (m/yr) 
% Ts = surface temperature (K)
% basal_slope & surface_slope = basal and surface slope effect (m/yr)
%      | ----- calculate by the slope (gradient MATLAB function) * ice velocity in X-Y direction
melt_range = -20:0.01:20; % range of basal melt rate (m/yr) 

% ---- 2. Simulating direct relationship between basal met rate and attenuation rate ----
% Initialization:
tmd.Na_SD = nan([length(melt_range) length(thickness)]); tmd.Na_TD = nan([length(melt_range) length(thickness)]); 
tmd.pts.Na_SD = nan([2 length(thickness)]); tmd.pts.Na_TD = nan([2 length(thickness)]); 
tmd.pts.m_SD = nan([2 length(thickness)]); tmd.pts.m_TD = nan([2 length(thickness)]); 

totalPoints = length(thickness);
h = waitbar(0, 'Simulating basal melt and attenuation relationships...');   % Initialize the waitbar
for n = 1:10%totalPoints
    if ~isnan(thickness(n)) && thickness(n)>100 && ~isnan(velocity(n)) && ...
            ~isnan(dTdl(n)) && ~isnan(smb(n)) && ~isnan(Ts(n)) && ~isnan(basal_slope(n)) && ~isnan(surface_slope(n)) 
        [tmd_Na_SD, tmd_Na_TD] = atten2melt(thickness(n), velocity(n), dTdl(n), smb(n), Ts(n), basal_slope(n), surface_slope(n), melt_range);
        tmd.Na_SD(:, n) = tmd_Na_SD';   % profiles with ice chemsitry from Siple Dome
        tmd.Na_TD(:, n) = tmd_Na_TD';   % profiles with ice chemsitry from Taylor Dome

        % ---- 3. Identifying boundary/inflection point (out-of-bound) ----
        [tmd.pts.Na_SD(:, n), tmd.pts.m_SD(:, n),] = off_bound_detection(tmd_Na_SD, melt_range);
        [tmd.pts.Na_TD(:, n), tmd.pts.m_TD(:, n),] = off_bound_detection(tmd_Na_TD, melt_range);
    end
    waitbar(n / totalPoints, h, sprintf('Progress: %.1f%%', n / totalPoints * 100));   % Update the progress bar and text
end
close(h);   % Close the progress bar when done

%% Calculate basal melt rate from radar-observed attenuation rates

% ---- 1. Loading attenuation rates ----
% For example:
load(dataFile, 'Attenuation', 'Attenuation_unc');   % results of attenuation_calculation.m

% ---- 2. Filtering/smoothing attenuation rate and uncertainty ----
Attenuation(Attenuation < 0) = NaN; Attenuation_unc(Attenuation < 0) = NaN;
Na = smoothdata(Attenuation, 250, 'gaussian', 'omitnan'); NaUnc = smoothdata(Attenuation_unc, 250, 'gaussian', 'omitnan');

% ---- 3. Interpolating basal melt rates ----
unc_SD = nan(size(Attenuation)); mean_SD = nan(size(Attenuation)); % Siple Dome chemistry
unc_TD = nan(size(Attenuation)); mean_TD = nan(size(Attenuation)); % Taylor Dome chemistry

for nn = 1:numel(Attenuation)
    valid = find(isfinite(tmd.Na_SD(:, nn)));
    [~,idx] = unique(tmd.Na_SD(:, nn));
    if ~isempty(valid) && length(idx) == length(melt_range)
        rand_Na = Na(nn) - NaUnc(nn) + NaUnc(nn) * 2 .* rand(100, 1); % generate 100 random values
        % Siple Dome:
        SD = interp1(tmd.Na_SD(valid, nn), melt_range(valid), rand_Na); 
        mean_SD(nn) = nanmean(SD); unc_SD(nn) = nanstd(SD)/sqrt(length(rand_Na)); 
        % Taylor Dome: 
        TD = interp1(tmd.Na_TD(valid, nn), melt_range(valid), rand_Na); 
        mean_TD(nn) = nanmean(TD); unc_TD(nn) = nanstd(TD)/sqrt(length(rand_Na)); 
    end
end
MeltRate.mu.SD = mean_SD; MeltRate.mu.TD = mean_TD;
MeltRate.unc.SD = unc_SD; MeltRate.unc.TD = unc_TD; 

% ---- 4. Filtering basal melt rates (if you want to) ----
MeltRate.mu.SD(MeltRate.mu.SD < tmd.pts.m_SD(1, :) | MeltRate.mu.SD > tmd.pts.m_SD(2, :)) = NaN;
MeltRate.mu.TD(MeltRate.mu.TD < tmd.pts.m_TD(1, :)  | MeltRate.mu.TD > tmd.pts.m_TD(2, :)) = NaN;
MeltRate.unc.SD(~isfinite(MeltRate.mu.SD)) = NaN;
MeltRate.unc.TD(~isfinite(MeltRate.mu.TD)) = NaN;

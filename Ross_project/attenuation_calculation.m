% ==== Author: Kiera Tran ====
% This script is used for calculating ice-penetrating radar attenuation
clear all; clc;

% add general directory to all existing functions in Ross_project folder
addpath(genpath('/Users/kieratran/Desktop/Research/Projects/ris/rosetta/Ross_project'))  

%% Load data:
% In this section, you can load radargram add any directory or folder and load data
% Best way to test run the code is to use Operation IceBridge data: https://data.cresis.ku.edu/data/rds/
%     |---- Note: use CSARP_standard both (L1B and L2B data). For example:
%           |---- echogram_fn = 'E:\rds\2012_Greenland_P3\CSARP_standard\20120514_01\Data_20120514_01_014.mat';
%           |---- layer_fn = 'E:\rds\2012_Greenland_P3\CSARP_layerData\20120514_01\Data_20120514_01_014.mat';
path = '/Users/kieratran/Desktop/Research/Projects/wais/2026/data/CreSIS/Venable/2016_DC8/';
dataFile = 'Data_20161104_05_033.mat';
echogram_fn = [path 'L1B/' dataFile];
layer_fn = [path 'L2B/' dataFile];

load(echogram_fn); load(layer_fn); % load all variables

% Neccessary variables: 'Data', 'surfLayer', 'bedLayer', 'Time', 'Latitude', 'Longitude', 'Elevation', 'surfTime', 'bedTime'
% Surface and Bed layers can either be in elevation, index, or time
% In this case, 'surfLayer' and 'bedLayer' are indices of surface and base layers
%     |---- Check units and conversions
%         |---- depth = -0.5 * (time - surface_time) * 1.68e8;   %  Sound velocity through ice is 1.68e8 m/s

surfTime = Surface; 
bedTime = layerData{2}.value{2}.data;
if ~exist("surfLayer")
    surfLayer = nan(size(Latitude));
    for n = 1:numel(Latitude)
        [~, surfLayer(n)] = nanmin(abs(Time - surfTime(n)));
    end
end

if ~exist("bedLayer")
    bedLayer = nan(size(Latitude));
    for n = 1:numel(Latitude)
        [~, bedLayer(n)] = nanmin(abs(Time - bedTime(n)));
    end
end


%% Find inflection points:
% This section identify (1) the boundary between firn and ice, 
% and (2) bed-echo length to neglect scattering power signals

warning off;

% ---- 1. Initialization ----
inLayer = nan(size(Latitude));   % indices of firn depth
inLayer2 = nan(size(Latitude));   % indices of bed-echo return
firn_confidence = nan(size(Latitude));

% ---- 2. Setting limits for each layers (you can tune this) ----
firn_prior = 70;   echo_prior = 50;   % expected depths (m)
firn_range = 30;   echo_range = 50;   % search ± range (m)
disp('>>>> Getting layers ...')
for i = 1:numel(Latitude)
    if isnan(surfLayer(i)) || isnan(bedLayer(i)) || bedLayer(i) <= surfLayer(i) + 10
        continue   % all outputs already NaN from pre-allocation
    end
    [power, depth] = geo_correction(Data(:,i), Time, surfTime(i), Elevation(i)); 

% ---- 3.  Dectecting layers' boundaries ----
    surf_s = surfLayer(i);
    bed_s = bedLayer(i);
    depth_col = depth(surf_s:bed_s);
    power_col = power(surf_s:bed_s);
    n_samp = numel(depth_col);
    dz = abs(mean(diff(depth_col)));
    % --- Firn/ice: search from surface downward ---
    firn_lo = max(2, round((firn_prior - firn_range) / dz));
    firn_hi = min(n_samp-1, round((firn_prior + firn_range) / dz));
    if firn_hi > firn_lo + 2
        [bp, firn_confidence(i)] = piecewise_fit(depth_col, power_col, firn_lo, firn_hi, "firn");
        inLayer(i) = surf_s + bp - 1;
    else
        warning('nn=%d: firn window too narrow (dz=%.2fm)', i, dz);
    end
    % --- Ice/saline: search from bed upward ---
    echo_lo = max(2, n_samp - round((echo_prior + echo_range) / dz));
    echo_hi = min(n_samp-1, n_samp - round((echo_prior - echo_range) / dz));
    if echo_hi > echo_lo + 2
        % Flip segment so "surface side" is always segment 1 in the fit
        d_echo = depth_col(echo_lo:echo_hi);
        p_echo = power_col(echo_lo:echo_hi);
        [bp, ~] = piecewise_fit(d_echo, p_echo, 2, numel(d_echo) - 1, "saline");
        % Map back — bp is relative to saline segment
        inLayer2(i) = surf_s + echo_lo + bp - 2;
    else
        warning('nn=%d: saline window too narrow (dz=%.2fm)', i, dz);
    end
end

% If you don't have 'smoothn' function install, download it here: 
%               https://www.mathworks.com/matlabcentral/fileexchange/25634-smoothn

% ---- 4. Smoothing layers ----
nanidx = isnan(inLayer); inLayer = smoothn(inLayer, 600); 
inLayer(nanidx) = NaN; inLayer = round(inLayer);
nanidx = isnan(inLayer2); inLayer2 = smoothn(inLayer2, 600); 
inLayer2(nanidx) = NaN; inLayer2 = round(inLayer2);

% ---- 5. Converting from indices to elevation (m) ----
firnElevation = nan(size(Latitude));
echoElevation = nan(size(Latitude));
surfElevation = nan(size(Latitude));
bedElevation = nan(size(Latitude));
for i = 1:numel(Latitude)
    if ~isnan(surfTime(i)) && ~isnan(inLayer(i)) && ~isnan(inLayer2(i))
        [power, depth] = geo_correction(Data(:,i), Time, surfTime(i), Elevation(i));
        firnElevation(i) = depth(inLayer(i));   % firn layer depth (m)
        echoElevation(i) = depth(inLayer2(i));   % bed-echo layer depth (m)
        surfElevation(i) = depth(surfLayer(i));   % surface layer depth (m)
        bedElevation(i) = depth(bedLayer(i));    % bed layer depth (m)
    end
end

% ---- 6. Calculating ice thickness (m) ----
IceThick = abs(bedElevation - surfElevation);

%% Calculate attenuation rate:
% This section calculates depth-resolved multi-reflector attenuation rates

% ---- 1. Initialization ----
AttenRate = nan(size(Latitude));   % linear fitting attenuation rate
AttenRate_unc = nan(size(Latitude));   % linear fitting uncertainty
Na = nan(size(Latitude));   % piecewise fitting attenuation rate
Na_unc = nan(size(Latitude));    % piecewise fitting uncertainty

disp('>>> Calculating attenuation rates ...')
for i = 1:numel(Latitude)
    %disp(i)
    if ~isnan(inLayer(i)) && ~isnan(inLayer2(i)) && inLayer2(i) > inLayer(i)   % surface and bed picks exist
        [power, depth] = geo_correction(Data(:,i), Time, surfTime(i), Elevation(i));   % geometrically correct returned power

% ---- 2. Looping through each segment length and calculate attenuation ----
        num = 1;
        for ii = 20:5:50 % segment length is based on percentage of total ice thickness (20% to 50% of ice thickness)
            % --- Piecewise attenuation rate
            atten_rate = atten_calc(power, depth, inLayer(i), inLayer2(i), 0, ii);
            atten_depth = depth_ave(atten_rate, depth, inLayer(i):inLayer2(i));
            Na(num, i) = atten_depth(inLayer2(i));
            % --- Piecewise uncertainty
            atten_rate = piecewise_unc(power, depth, inLayer(i), inLayer2(i), 0, ii);
            atten_depth = depth_ave(atten_rate, depth, inLayer(i):inLayer2(i));
            Na_unc(num, i) = atten_depth(inLayer2(i));
            num = num+1;
        end
        % --- Segment length = 100% ice thickness & uncertainty
        p = polyfit(depth(inLayer(i):inLayer2(i)), power(inLayer(i):inLayer2(i)),1)/2*1000;   % One-way linear fitting attenuation rate (dB/km)
        AttenRate(i) = p(1,1);
        AttenRate_unc(i) = slopeSE(depth(inLayer(i):inLayer2(i)), power(inLayer(i):inLayer2(i)));   % uncertainty
        % --- Combine all segment lengths
        Na(num, i) = AttenRate(i);
        Na_unc(num, i) = AttenRate_unc(i);
    end
end

% ---- 3. Averaging depth-resolved attenuation rates with different piecewise regressions
 Attenuation = nanmean(Na);
 Attenuation_unc = nanmean(Na_unc);


%% Sanity check and visualization:
% This will load A-scope (radargram), 
% and you can choose specific Z-scope (power profile) to view for validation
plot_radargram
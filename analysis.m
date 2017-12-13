% -------------------------------------------------------------------------
% This program analyses the results generated in the main program
% 
% USAGE
% ---------
% load the 'RotorArray' and 'RUN_TYPE' properties of the workspace from 'main.m'
% WORKSPACES is the array of workspaces that have to be compared
% RUN_TYPE will be dispalyed as legend
% -------------------------------------------------------------------------

%% setting up
clear;
close all;
clc;

RUN_DESC = string('Correction');
SAVE_FIG = 0;

disp(sprintf('Analysing results...')); %#ok<*DSPS>

%% load .mat files
disp(sprintf('Loading workspaces'));
% Betz3B          = load('results/Betz_0929_2012.mat', 'RotorArray', 'RUN_TYPE');
% Betz1B          = load('results/Betz1B_0924_2231.mat', 'RotorArray', 'RUN_TYPE');
% Betz2B          = load('results/Betz2B_0924_2234.mat', 'RotorArray', 'RUN_TYPE');
% Betz4B          = load('results/Betz4B_0924_2211.mat', 'RotorArray', 'RUN_TYPE');
% Betz5B          = load('results/Betz5B_0924_2227.mat', 'RotorArray', 'RUN_TYPE');
 BetzNoCorr      = load('results/BetzP0_1011_1639.mat', 'RotorArray', 'RUN_TYPE');
% BetzSimple      = load('results/BetzSimple_0929_2012.mat', 'RotorArray', 'RUN_TYPE');
% BetzPolyfit     = load('results/BetzPolyfit_0920_1556.mat', 'RotorArray', 'RUN_TYPE');
% BetzBezfit      = load('results/BetzBezierfit_0927_1802.mat', 'RotorArray', 'RUN_TYPE');
% Emperical12B    = load('results/Emperical12B_0924_2246.mat', 'RotorArray', 'RUN_TYPE');
% Emperical13B    = load('results/Emperical13B_0924_2247.mat', 'RotorArray', 'RUN_TYPE');
% Emperical21B    = load('results/Emperical2-1B_0925_2331.mat', 'RotorArray', 'RUN_TYPE');
% Emperical22B    = load('results/Emperical22B_0924_2247.mat', 'RotorArray', 'RUN_TYPE');
% Emperical23B    = load('results/Emperical23B_0924_2247.mat', 'RotorArray', 'RUN_TYPE');
% Emperical24B    = load('results/Emperical2-4B_0925_2331.mat', 'RotorArray', 'RUN_TYPE');
% Emperical25B    = load('results/Emperical2-5B_0925_2331.mat', 'RotorArray', 'RUN_TYPE');
% PSOBezier       = load('results/PSOBezier_0929_1923.mat', 'RotorArray', 'RUN_TYPE');
% PSOPolynomial   = load('results/PSOPolynomial_0928_2300.mat', 'RotorArray', 'RUN_TYPE');
% DEOBezier       = load('results/DEOBezier_0929_1831.mat', 'RotorArray', 'RUN_TYPE');
% DEOPolynomial   = load('results/DEOPolynomial_0928_2157.mat', 'RotorArray', 'RUN_TYPE');
% POn_GOn         = load('results/PraOn-GlaOn_0925_1409.mat', 'RotorArray', 'RUN_TYPE');
% POn_GOff        = load('results/PraOn-GlaOff_0924_2319.mat', 'RotorArray', 'RUN_TYPE');
% POff_GOn        = load('results/PraOff-GlaOn_0924_2320.mat', 'RotorArray', 'RUN_TYPE');
% POff_GOff       = load('results/PraOff-GlaOff_0924_2319.mat', 'RotorArray', 'RUN_TYPE');
% BetzDU          = load('results/DU-95-W-180_0925_1521.mat', 'RotorArray', 'RUN_TYPE');
% BetzNACA        = load('results/NACA-63-415_0925_1959.mat', 'RotorArray', 'RUN_TYPE');
% BetzAH          = load('results/AH-93-W-257_0925_1959.mat', 'RotorArray', 'RUN_TYPE');
% Emp1DU          = load('results/Emp1_DU-95-W-180_0925_2011.mat', 'RotorArray', 'RUN_TYPE');
% Emp1NACA        = load('results/Emp1_NACA-63-415_0925_2011.mat', 'RotorArray', 'RUN_TYPE');
% Emp1AH          = load('results/Emp1_AH-93-W-257_0925_2010.mat', 'RotorArray', 'RUN_TYPE');
% Emp2AH          = load('results/Emperical2AH_0929_2122.mat', 'RotorArray', 'RUN_TYPE');
% BetzSmallHub    = load('results/BetzSmallHub_0929_2117.mat', 'RotorArray', 'RUN_TYPE');
% BetzVSmallHub   = load('results/BetzVSmallHub_0929_2116.mat', 'RotorArray', 'RUN_TYPE');
Pitch_1         = load('results/BetzP_1_1011_1639.mat', 'RotorArray', 'RUN_TYPE');
Pitch3         = load('results/BetzP3_1011_1639.mat', 'RotorArray', 'RUN_TYPE');
Pitch1          = load('results/BetzP1_1011_1639.mat', 'RotorArray', 'RUN_TYPE');
Pitch2          = load('results/BetzP2_1011_1639.mat', 'RotorArray', 'RUN_TYPE');
% BetzSmallAH     = load('results/BetzSmallAH_0929_2137.mat', 'RotorArray', 'RUN_TYPE');
% BetzVSmallAH    = load('results/BetzVSmallAH_0929_2137.mat', 'RotorArray', 'RUN_TYPE');
% BetzSmall2B     = load('results/BetzSmall2B_0929_2154.mat', 'RotorArray', 'RUN_TYPE');
% BetzVSmall2B    = load('results/BetzVSmall2B_0929_2154.mat', 'RotorArray', 'RUN_TYPE');
% BetzNoHub       = load('results/BetzNoHub_0929_2204.mat', 'RotorArray', 'RUN_TYPE');
disp(sprintf('Loading complete'));


%% Include the RUN_TYPE that you need to compare
Workspaces      = [Pitch_1 BetzNoCorr Pitch1 Pitch2 Pitch3];

%% print best Cp,max for each run type
for ws = Workspaces
    best_cp = max([ws.RotorArray(:).cP]);
    disp(sprintf('Run Type = %s, Cp,max = %.4f', ws.RUN_TYPE, best_cp));
end

%% CP-Lambda Curve
disp(sprintf('Generating Cp-Lambda curve'));
f = figure(1); 
%subplot(1,2,1)
for ws = Workspaces
    [best_cp, best_i]   = max([ws.RotorArray(:).cP]);
    BestRotor           = ws.RotorArray(best_i);
    peak                = mat2str([BestRotor.lambda round(BestRotor.cP,4)]);
    
    plot([ws.RotorArray(:).lambda], [ws.RotorArray(:).cP], 'DisplayName', char(ws.RUN_TYPE));
        xlabel('Tip Speed Ratio [-]');
        ylabel('Power Coefficient [-]');
        hold on;
        
    %text(BestRotor.lambda,  BestRotor.cP,  peak);  % peak of the curve
       
end     
legend('show');
hold off;
if(SAVE_FIG) 
    saveas(f, sprintf('plots/cp_lambda_%s_%s.png', RUN_DESC, datestr(now,'mmdd_HHMM'))); 
end
        
%% Blade Design
disp(sprintf('Generating design profile'));
f = figure(2);

% chord
subplot(1,2,1);
for ws = Workspaces
    [best_cp, best_i] = max([ws.RotorArray(:).cP]);
    BestRotor        = ws.RotorArray(best_i);
    plot([BestRotor.Annuli.mu], [BestRotor.Annuli.c], 'DisplayName', char(ws.RUN_TYPE));
        xlabel('r/R [-]');
        ylabel('Chord Length [m]');
        hold on;
end 
legend('show');
hold off;

% twist
subplot(1,2,2);
for ws = Workspaces
    [best_cp, best_i]   = max([ws.RotorArray(:).cP]);
    BestRotor           = ws.RotorArray(best_i);
    label               = string(ws.RUN_TYPE);
    plot([BestRotor.Annuli.mu], [BestRotor.Annuli.twist], 'DisplayName', char(ws.RUN_TYPE));
        xlabel('r/R [-]');
        ylabel('Twist Angle [deg]');
        hold on;
end 
legend('show');
hold off;

if(SAVE_FIG) 
    saveas(f, sprintf('plots/design_%s_%s.png', RUN_DESC, datestr(now,'mmdd_HHMM')));
end



%% Spanwise distribution
fields = {...
%           'c',      'Chord Length [m]'; ...
%           'twist',  'Twist Angle [degree]'; ...
%           'alpha',  'Angle of Attack [degree]'; ...
%           'phi',    'Inflow Angle [degree]'; ...
%           'f',      'Tip-Root Correction [-]'; ...
%           'aA',     'Axial Induction Factor [-]'; ...
%           'aT',     'Tangential Induction Factor [-]'; ...
%           'cX',     'Axial Force Coefficient [-]'; ...
%           'cY',     'Azimuthal Force Coefficient [-]'; ...
%           'cQ',     'Torque Coefficient [-]';...
%           'cT',     'Thrust Coefficient [-]';...
            'cP',     'Power Coefficient [-]';...
          };
size_fields = size(fields);

for i = 1:size_fields(1,1)
    f = figure(i+2);
    
    for ws = Workspaces
        [best_cp, best_i]   = max([ws.RotorArray(:).cP]);
        %best_i              = 8;
        BestRotor           = ws.RotorArray(best_i);
        
        plot([BestRotor.Annuli(:).mu], [BestRotor.Annuli(:).(fields{i,1})], 'DisplayName', char(ws.RUN_TYPE));
            xlabel('r/R [-]');
            ylabel(fields{i,2});
            hold on;
    end 
    
    legend('show');
    hold off;
    
    if(SAVE_FIG) 
        saveas(f, sprintf('plots/%s_%s_%s.png', RUN_DESC, fields{i,1}, datestr(now,'mmdd_HHMM')));
    end
end
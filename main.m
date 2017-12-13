% -------------------------------------------------------------------------
% TU DELFT - System Integration Project 
% 
% Created by    - ttanmay@tudelft.nl
% Supervised by - M.B.Zaayer@tudelft.nl
% Created on    - 13th August, 2017
%
% The objective of this code is to :
% run BEM code under various options
%
% ------------------------------------------
% Acceptable values for BEM design 'profile'
% ------------------------------------------
% 'BetzSimple'      design parameters determined using Betz limit of aA=1/3 
%                   without including wake rotation
% 'Betz'            design parameters determined using Betz limit of aA=1/3 
%                   including wake rotation
% 'Polynomial'      design parameters as polynomial of degree 'N_POLYFIT'
% 'Bezier'          design parameters as bezier with 3 control points 
%
% ------------------------------------------
% Acceptable values for RUN_TYPE
% ------------------------------------------
% 'BetzSimple'      design parameters determined using Betz limit of aA=1/3 
%                   without including wake rotation
% 'Betz'            design parameters determined using Betz limit of aA=1/3 
%                   including wake rotation
% 'BetzNoCorr'      run 'Betz' without Prandtl & Glauert correction
% 'BetzPolyfit'     polyfit of degree 'N_POLYFIT' for profile obtained from
%                   'Betz'
% 'BetzBezierfit'   bezierfit with 3 control points for profile obtained from
%                   'Betz'
% 'PSOPolynomial'   PSO with design parameters as polynomial of degree 'N_POLYFIT'
% 'PSOBezier'       PSO with design parameters as bezier with 3 control points 
% 'DEOPolynomial'   DEO with design parameters as polynomial of degree 'N_POLYFIT'
% 'DEOBezier'       DEO with design parameters as bezier with 3 control points
% 'Emperical1'      Emperical relationship of cp and lambda. Wind Energy Conversion Systems: Freris (1990)
% 'Emperical2'      Emperical relationship of cp and lambda. Wind Enrgy Explained: Manwell (2010)
% -------------------------------------------------------------------------

%% Setting Up
close all;
clear;
clc;


%% Model Parameters
lambda_array    = 5:0.5:12;                 % sample space for tip speed ratios [-]
N_POLYFIT       = 3;                        % degree of polynomial function
RUN_TYPE        = string('Betz');      % refer RUN_TYPE in code header
RUN_DESC        = string('BetzP3');      % description that will appear in legends

n_lambda        = numel(lambda_array);
t_main          = tic;

disp(sprintf('Started running %s...', RUN_TYPE)); %#ok<*DSPS>

%% Variables to analyse
RotorArray(n_lambda, 1)    = RotorClass();  % to store data for rotor optimised for a given lambda

%% Evaluating various rotor designs
 for i = 1:n_lambda       
        lambda = lambda_array(i);
        
        disp(sprintf('Evaluating TSR = %.1f', lambda));
        
        %% Betz's analytical optimization with wake rotation
        % BEMRotor(lambda, pitch, Xc, Xt, profile, isPrandtl, isGlauert)
        Rotor = BEMRotor(lambda, 3, [], [], string('Betz'), 1, 1);
        
        %% Betz's analytical optimization without wake rotation
        if(RUN_TYPE == 'BetzSimple')
            Rotor = BEMRotor(lambda, 0, [], [], string('BetzSimple'), 1, 1); 
        end
        
        %% Polyfit to Betz's optimal
        if(RUN_TYPE == 'BetzPolyfit')
            XC_INIT = polyfit([Rotor.Annuli(:).mu], [Rotor.Annuli(:).c], N_POLYFIT);
            XT_INIT = polyfit([Rotor.Annuli(:).mu], [Rotor.Annuli(:).twist], N_POLYFIT);
            Rotor = BEMRotor(lambda, 0, XC_INIT, XT_INIT, string('Polynomial'), 1, 1); 
        end
        
        %% Bezierfit to Betz's optimal
        if(RUN_TYPE == 'BetzBezierfit')
            mu = [Rotor.Annuli.mu];
            c = [Rotor.Annuli.c];
            t = [Rotor.Annuli.twist];
            XC_INIT = [c(1) mu(25) c(25) c(50)];
            XT_INIT = [t(1) mu(25) t(25) t(50)];
            
            Rotor = BEMRotor(lambda, 0, XC_INIT, XT_INIT, string('Bezier'), 1, 1); 
        end
        
        %% 'Betz' without Prandtl & Glauert correction
        if(RUN_TYPE == 'BetzNoCorr')
            Rotor = BEMRotor(lambda, 0, [], [], string('Betz'), 0, 0);
        end
        
        %% emperical relationship. Freris (1990)
        if(RUN_TYPE == 'Emperical1')
            CL_CD_MAX   = max(Rotor.airfoil(:,2)./Rotor.airfoil(:,3));
            Rotor.cP = (16/27)*((1-(0.416/(Rotor.n_blades*lambda)))^2)*(exp(-0.35/(lambda^1.29))-lambda/CL_CD_MAX);
        end
        
        %% emperical relationship. Manwell (2010)
        if(RUN_TYPE == 'Emperical2')
            CL_CD_MAX   = max(Rotor.airfoil(:,2)./Rotor.airfoil(:,3));
            Rotor.cP = (16*lambda/27)*(lambda + (1.32 + (lambda-8)^2/(400))/(Rotor.n_blades^(2/3)))^-1 - (0.57*lambda^2)/((CL_CD_MAX)*(lambda+1/6));
        end 
        
        %% Particle Swarm Optimization with polynomial parameterization        
        if(RUN_TYPE == 'PSOPolynomial')
            XC_INIT = polyfit([Rotor.Annuli(:).mu], [Rotor.Annuli(:).c], N_POLYFIT);
            XT_INIT = polyfit([Rotor.Annuli(:).mu], [Rotor.Annuli(:).twist], N_POLYFIT);
            
            BestDesign = ParticleSwarm(lambda, 0, XC_INIT, XT_INIT, string('Polynomial'), 1, 1);
            Rotor = BEMRotor(lambda, 0, BestDesign.Xc_best, BestDesign.Xt_best, string('Polynomial'), 1, 1);
        end
        
        %% Particle Swarm Optimization with bezier parameterization     
        if(RUN_TYPE == 'PSOBezier')
            % initialize Xc & Xt for Bezier 
            mu = [Rotor.Annuli.mu];
            c = [Rotor.Annuli.c];
            t = [Rotor.Annuli.twist];
            XC_INIT = [c(1) mu(25) c(25) c(50)];
            XT_INIT = [t(1) mu(25) t(25) t(50)];
            
            BestDesign = ParticleSwarm(lambda, 0, XC_INIT, XT_INIT, string('Bezier'), 1, 1);
            Rotor = BEMRotor(lambda, 0, BestDesign.Xc_best, BestDesign.Xt_best, string('Bezier'), 1, 1);
        end   
        
        %% Differential Evolution Optimization with polynomial parameterization     
        if(RUN_TYPE == 'DEOPolynomial')
            XC_INIT = polyfit([Rotor.Annuli(:).mu], [Rotor.Annuli(:).c], N_POLYFIT);
            XT_INIT = polyfit([Rotor.Annuli(:).mu], [Rotor.Annuli(:).twist], N_POLYFIT);
            
            BestDesign = DifferentialEvolution(lambda, 0, XC_INIT, XT_INIT, string('Polynomial'), 1, 1);
            Rotor = BEMRotor(lambda, 0, BestDesign.Xc, BestDesign.Xt, string('Polynomial'), 1, 1);
        end
        
        %% Differential Evolution Optimization with bezier parameterization    
        if(RUN_TYPE == 'DEOBezier')
            % initialize Xc & Xt for Bezier 
            mu = [Rotor.Annuli.mu];
            c = [Rotor.Annuli.c];
            t = [Rotor.Annuli.twist];
            XC_INIT = [c(1) mu(25) c(25) c(50)];
            XT_INIT = [t(1) mu(25) t(25) t(50)];
            
            BestDesign = DifferentialEvolution(lambda, 0, XC_INIT, XT_INIT, string('Bezier'), 1, 1);
            Rotor = BEMRotor(lambda, 0, BestDesign.Xc, BestDesign.Xt, string('Bezier'), 1, 1);
        end  
        
        RotorArray(i) = Rotor;
 end % end of lambda loop

[~, i] = max([RotorArray(:).cP]);
BestRotor = RotorArray(i);

%% Saving
t_elapsed   = toc(t_main);
RUN_TYPE    = RUN_DESC;
save(sprintf('results/%s_%s.mat', RUN_TYPE, datestr(now,'mmdd_HHMM')));
disp(sprintf('Completed in %.1f min.', t_elapsed/60));


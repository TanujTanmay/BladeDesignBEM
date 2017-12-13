% -------------------------------------------------------------------------
% The objective of this function is to :
% perform the BEM modeling over the entire rotor

% INPUTS
% --------------
% lambda        design tip speed ratio of the rotor     [-]
% pitch         design default blade pitch              [deg]
% Xc            chord length coefficients               [aN ... a2 a1 a0]
% Xt            twist angle coefficients                [aN ... a2 a1 a0]
% profile       profile for chord & twist distribution  [string]
% isPrandtl     is Prandtl correction active            [0/1]
% isGlauert     is Glauert correction active            [0/1]
%
% OUTPUT
% ---------------
% Rotor         properties of the rotor                 [RotorClass]
%
% CHANGE LOGS
% ---------------
% 07 Sep 2017   included polynomial coefficients for chord & twist distribution in the Rotor properties 
% 09 Sep 2017   handling bezier and polynomial functions
% 12 Sep 2017   computational time is now recorded
% 19 Sep 2017   parametrized profile and corrections
% -------------------------------------------------------------------------

function Rotor = BEMRotor(lambda, pitch, Xc, Xt, profile, isPrandtl, isGlauert)

    t_rotor = tic;
    Rotor = RotorClass(lambda, pitch, Xc, Xt, profile, isPrandtl, isGlauert);  
    
    %% determine optimal airfoil properties
    [~, i_opt]   = max(Rotor.airfoil(:,2)./Rotor.airfoil(:,3));
    ALPHA_OPT    = Rotor.airfoil(i_opt,1);
    CL_OPT       = Rotor.airfoil(i_opt,2);
    
    %% compute bezier functions only when required
    if(profile == 'Bezier')
        BezierChord  = bezier(Xc);
        BezierTwist  = bezier(Xt);
    end
    
    %% split rotor into annuli
    r_array = cosspace(Rotor.root, Rotor.radius, Rotor.n_annuli+1); 
    Rotor.Annuli(Rotor.n_annuli,1) = AnnulusClass();
    
    %% evaluate each annulus
    for i = 1:Rotor.n_annuli
        
        %% initialize the annulus geometry
        Annulus          = AnnulusClass();
        Annulus.r1       = r_array(i);
        Annulus.r2       = r_array(i+1);
        Annulus.r        = 0.5*(Annulus.r1 + Annulus.r2);
        Annulus.mu       = Annulus.r/Rotor.radius;
        Annulus.dr       = Annulus.r2 - Annulus.r1;
        Annulus.area     = 2*pi*Annulus.r*Annulus.dr;
        Annulus.lambda_r = lambda*Annulus.mu;
        
        %% initialize twist and chord distribution
        if(profile == 'Bezier')
            Annulus.c      = interp1(BezierChord(:,1), BezierChord(:,2), Annulus.mu);
            Annulus.twist  = interp1(BezierTwist(:,1), BezierTwist(:,2), Annulus.mu);
        elseif(profile == 'Polynomial')
            Annulus.c      = polyval(Xc, Annulus.mu);
            Annulus.twist  = polyval(Xt, Annulus.mu);
        elseif(profile == 'Betz')
            % with wake rotation
            phi_opt        = (2/3) * atand(1/Annulus.lambda_r);
            Annulus.c      = (8*pi*Annulus.r*(1-cosd(phi_opt)))/(Rotor.n_blades*CL_OPT);
            Annulus.twist  = phi_opt - ALPHA_OPT - pitch; 
        elseif(profile == 'BetzSimple')
            % without wake rotation
            phi_opt        = atand(2/(3*Annulus.lambda_r));
            Annulus.c      = (8*pi*Rotor.radius*sind(phi_opt))/(3*Rotor.n_blades*CL_OPT*lambda);
            Annulus.twist  = phi_opt - ALPHA_OPT - pitch;
        else
            Annulus.c      = 1 + 3*(1-Annulus.mu);
            Annulus.twist  = 14*(1-Annulus.mu);
        end 

        %% run BEM over the annulus   
        Annulus            = BEMAnnulus(Annulus, Rotor, isPrandtl, isGlauert);
        Annulus.sigma_r    = (Rotor.n_blades*Annulus.c)/(2*pi*Annulus.r);
        Rotor.Annuli(i)    = Annulus;    

        %% compute the values at rotor level
        Rotor.cT = Rotor.cT + Annulus.cT*Annulus.area/Rotor.area;
        Rotor.cQ = Rotor.cQ + Annulus.cQ*Annulus.area/Rotor.area;
        Rotor.cP = Rotor.cP + Annulus.cP*Annulus.area/Rotor.area;
    end % end of BEM loop

    Rotor.t      = toc(t_rotor);
    
end
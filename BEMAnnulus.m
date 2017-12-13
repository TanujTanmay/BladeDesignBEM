% -------------------------------------------------------------------------
% The objective of this function is to :
% perform the BEM modeling over an input Annulus of defined rotor properties

% INPUTS
% --------------
% Annulus       geometry of the annulus                   [AnnulusClass]
% Rotor         properties of the rotor                   [RotorClass]
% isPrandtl     is Prandtl correction active              [0/1]
% isGlauert     is Glauert correction active              [0/1]
%
% OUTPUT
% ---------------
% Annulus       properties of the annulus                 [AnnulusClass]
%
% CHANGE LOGS
% ---------------
% 14 Aug 2017   removed PITCH from twist design
% 27 Aug 2017   included blade design with wake rotation
% 07 Sep 2017   parametrised chord and twist distribution
% 09 Sep 2017   handling bezier & polynomial functions for chord & twist
% 12 Sep 2017   computational time is now recorded
% -------------------------------------------------------------------------

function Annulus = BEMAnnulus(Annulus, Rotor, isPrandtl, isGlauert)

    t_annulus = tic;
    
    %% Global Variables
    DENSITY     = 1.225;            % air density                   [kg/m3]
    U_INFINITY  = 10;               % free wind speed               [m/s]
    CT1         = 1.816;            % thrust coefficient at aA=1    [-]    
    CT2         = 2*sqrt(CT1)-CT1;  % 2SQRT(CT1)-CT1                [-]

    %% Declare Iteration Parameters
    N_ITR       = 100;              % maximum number of iterations  [-]
    ITR_TOL     = 0.0001;           % tolerance for convergance     [-]
    ITR_RELAX   = 0.75;             % relax factor for iteration    [-]

    %% Retreive Rotor Parameters
    RADIUS      = Rotor.radius;
    ROOT        = Rotor.root/Rotor.radius;    
    N_BLADES    = Rotor.n_blades;
    PITCH       = Rotor.pitch;
    LAMBDA      = Rotor.lambda;
    AIRFOIL     = Rotor.airfoil;
    
    %% Retreive Annulus Parameters
    r           = Annulus.r;
    mu          = Annulus.mu;
    dr          = Annulus.dr;
    area        = Annulus.area;
    lambda_r    = Annulus.lambda_r;
    c           = Annulus.c;
    twist       = Annulus.twist;

    %% estimates of axial and tangential induction factors
    aA = 0.3;
    aT = 0;       
    
    %% iterate until axial induction factor converge
    for i = 1:N_ITR
        
        phi     = atand((1-aA)/(lambda_r*(1+aT)));
        alpha   = phi - twist - PITCH;        
              
        cL      = interp1(AIRFOIL(:,1),AIRFOIL(:,2),alpha);
        cD      = interp1(AIRFOIL(:,1),AIRFOIL(:,3),alpha);  
        
        % bound cL and cD
        if(cL < 1e-6 || isnan(cL))
            cL = 1e-6;
            cD = 0;
        end
        
        gamma   = 0.5*cL*U_INFINITY*c;
        w       = sqrt(((U_INFINITY*(1-aA))^2) + ((lambda_r*U_INFINITY*(1+aT))^2));
        lift    = 0.5*c*DENSITY*(w^2)*cL;
        drag    = 0.5*c*DENSITY*(w^2)*cD;
        
        fX      = (lift*cosd(phi))+(drag*sind(phi));
        fY      = (lift*sind(phi))-(drag*cosd(phi));
        cX      = fX/(0.5*DENSITY*(U_INFINITY^2)*RADIUS);
        cY      = fY/(0.5*DENSITY*(U_INFINITY^2)*RADIUS);
        cT      = (fX*N_BLADES*dr)/(0.5*DENSITY*(U_INFINITY^2)*area);
        cQ      = (fY*N_BLADES*dr*mu)/(0.5*DENSITY*(U_INFINITY^2)*area);        
        cP      = cQ*LAMBDA;      
        
        %% Prandtl correction for tip and root losses        
        if(isPrandtl)
            f_tip   = (2/pi)*acos(exp(-(N_BLADES/2)*((1-mu)/mu)*sqrt(1+((LAMBDA*mu)/(1-aA))^2)));
            f_root  = (2/pi)*acos(exp(-(N_BLADES/2)*((mu-ROOT)/mu)*sqrt(1+((LAMBDA*mu)/(1-aA))^2)));
            f       = f_tip*f_root;
            f(f<0.0001)=0.0001;    
        else
            f_tip   = 1;
            f_root  = 1;
            f       = f_tip*f_root;
        end

        %% Glauert correction for heavily loaded rotor
        if(isGlauert)
            if (cT < CT2)
                aA_new = 0.5 - 0.5*sqrt(1 - cT);
            else
                aA_new = 1 + 0.25*(cT - CT1)/(sqrt(CT1) - 1);
            end
        else
            if(cT <= 0.96)
                aA_new  = 0.5 - 0.5*sqrt(1 - cT);
            else
                aA_new  = 0.4;
            end
        end
        
        aA_new = aA_new/f;
        
        %% Bound the value of axial induction
        if (aA_new > 0.96)
            aA_new = 0.96;
            cT = CT1 - 4*(sqrt(CT1) - 1)*(1-aA_new);
        end
        
        aA = ITR_RELAX*aA + (1-ITR_RELAX)*aA_new;
        
        aT = (fY*N_BLADES)/(2 * 2*pi*r * DENSITY * (U_INFINITY^2) * (1-aA) * lambda_r);
        aT = aT/f;
        
        %% Bound the value of tangential induction
        if (aT > aA*(1-aA)/(lambda_r^2))
            aT = aA*(1-aA)/(lambda_r^2);
        end     
        
        %% Check convergence
        if(abs((aA_new - aA)) <= ITR_TOL)
            break;
        end        
         

    end % end of iteration
    
    %% bound chord and twist distribution
    % or exclude the annulus that does not converge
    if (c < 0 || c > 30 || abs(twist) > 30 || i == N_ITR)
        cT = 0;
        cQ = 0;
        cP = 0;
    end 
    
    %% Save the results
    Annulus.phi     = phi;
    Annulus.alpha   = alpha;
    Annulus.gamma   = gamma;
    Annulus.aA      = aA;
    Annulus.aT      = aT;
    Annulus.w       = w;
    Annulus.f_tip   = f_tip;
    Annulus.f_root  = f_root;
    Annulus.f       = f;
    Annulus.cL      = cL;
    Annulus.cD      = cD;
    Annulus.lift    = lift;
    Annulus.drag    = drag;
    Annulus.fX      = fX;
    Annulus.fY      = fY;
    Annulus.cX      = cX;
    Annulus.cY      = cY;
    Annulus.cT      = cT;
    Annulus.cQ      = cQ;
    Annulus.cP      = cP;
    Annulus.conv    = i;  
    Annulus.t       = toc(t_annulus);
        
end % end of BEM function


classdef RotorClass
     
    %% Define the properties of the Rotor   
    properties        
        n_blades    = 3                         % number of blades              [-]        
        radius      = 50                        % radius of the rotor           [m]
        root        = 10                        % root of the blade             [m]        
        area        = pi*50^2                   % rotor swept area              [m^2]
        airfoil     = xlsread('DU-95-W-180.xlsx')   % airfoil properties        [excel]
        n_annuli    = 50                        % number of annulus             [-]
        lambda                                  % tip speed ratio               [-]
        pitch                                   % blade pitch angle             [deg] 
        Xc                                      % chord length coefficients     [aN ... a2 a1 a0]
        Xt                                      % twist angle coefficients      [aN ... a2 a1 a0]
        profile                                 % design parameter profile      [string]
        isPrandtl                               % is Prandtl correction active  [0/1] 
        isGlauert                               % is Glauert correction active  [0/1]    
        Annuli      = AnnulusClass();           % array of annulus              [AnnulusClass]
        cT = 0                                  % thrust force coefficient      [-]            
        cQ = 0                                  % torque coefficient            [-]
        cP = 0                                  % power coefficient             [-]  
        t                                       % computational time            [s]        
    end     
    
    %% Constructor class
    methods        
        function Rotor = RotorClass(lambda, pitch, Xc, Xt, profile, isPrandtl, isGlauert)            
            if nargin == 0
                % skipping error handling
            else
                Rotor.lambda  = lambda;
                Rotor.pitch   = pitch;
                Rotor.Xc      = Xc;
                Rotor.Xt      = Xt;
                Rotor.profile = profile;
                Rotor.isPrandtl = isPrandtl;
                Rotor.isGlauert = isGlauert;                
            end
        end    
    end
    
end     % end of class  
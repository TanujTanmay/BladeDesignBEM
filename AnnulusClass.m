classdef AnnulusClass
     
    %% Define the properties of the Annulus   
    properties
        r1                  % starting radial position      [m]
        r2                  % ending radial position        [m]
        r                   % mean radial position          [m]
        mu                  % non-dimensional r             [-]
        dr                  % thickness of annuli           [m]
        area                % area of annuli                [m^2]
        c                   % chord length                  [m]
        sigma_r             % local solidity                [-]
        lambda_r            % local speed ratio             [-]    
        phi                 % inflow angle                  [degree]
        twist               % section twist angle           [degree]
        alpha               % angle of attack               [degree]
        gamma               % circulation                   [m^2/s]
        aA                  % axial induction factor        [-]        
        aT                  % tangential induction factor   [-]
        w                   % local apparent wind speed     [m/s]
        f_tip               % tip loss factor               [-]
        f_root              % root loss factor              [-]
        f                   % f_tip*f_root                  [-]
        cL                  % lift coefficient              [-]
        cD                  % drag coefficient              [-]
        lift                % lift force                    [N/m]
        drag                % drag force                    [N/m]
        fX                  % 2D axial force                [N/m]
        fY                  % 2D azimuthal force            [N/m]
        cX                  % non-dimensional fX            [-]
        cY                  % non-dimensional fY            [-]
        cT                  % thrust force coefficient      [-]
        cQ                  % torque coefficient            [-]
        cP                  % power coefficient             [-]
        conv                % number of iterations          [-]
        t                   % computational time            [s]
    end     
    
end     % end of class  



    
        
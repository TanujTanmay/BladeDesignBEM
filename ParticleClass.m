classdef ParticleClass
     
    %% Define the properties of Particles in the Swarm
    properties        
        Xc          % chord coefficients            [aN ... a0]
        Xc_best     % best chord coefficients       [aN ... a0]
        Xc_velocity % velocity of chord particle    [aN ... a0]        
        Xt          % twist coefficients            [aN ... a0]
        Xt_best     % best twist coefficients       [aN ... a0]
        Xt_velocity % velocity of twist particle    [aN ... a0]        
        cP          % cp of the particle            [-]  
        t           % optimization time             [s]
    end     
    
    %% Constructor class
    methods        
        function obj = ParticleClass(id, XC_INIT, XT_INIT) 
            
            if nargin == 0
                % skipping error handling
            else
                if(id == 1)
                    obj.Xc = XC_INIT;   % at least 1 particle should have optimal XC and XT
                    obj.Xt = XT_INIT;   
                else
                    obj.Xc = (diag(2 - 4*rand(numel(XC_INIT),1))*XC_INIT')';    % -2 to +2 times of each coefficient of XC_INIT
                    obj.Xt = (diag(2 - 4*rand(numel(XT_INIT),1))*XT_INIT')';    % -2 to +2 times of each coefficient of XT_INIT
                end       
                obj.Xc_best     = obj.Xc;
                obj.Xc_velocity = zeros(1, numel(XC_INIT));          
                obj.Xt_best     = obj.Xt;
                obj.Xt_velocity = zeros(1, numel(XT_INIT));    
            end
            
        end % function   
    end % methods
    
end     % end of class  
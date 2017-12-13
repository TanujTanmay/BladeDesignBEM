classdef IndividualClass
     
    %% Define the properties of an Individual in the Population
    properties
        Xc      % coefficients of chord distribution    [aN ... a0]
        Xt      % coefficients of twist distribution    [aN ... a0]
        cP      % value of cp at new position           [-]
        t       % optimization time                     [s]
    end     
    
    %% Constructor class
    methods        
        function obj = IndividualClass(id, XC_INIT, XT_INIT) 
            
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
            end
            
        end % function   
    end % methods
    
end     % end of class
% -------------------------------------------------------------------------
% The objective of this function is to :
% create a Bezier curve based on input hinged points
%
% Hinged Points - (0.2, X(0)); (X(1), X(2)); (X(3), X(4)); ...; (1, X(N))
%
% INPUTS
% --------------
% X            coefficient for hinged points                 [aN ... a2 a1 a0]
%
% OUTPUT
% ---------------
% Bezier       range of values over certain domain           [x f(x)]
%
% CHANGE LOGS
% ---------------
% 09 Sep 2017   created
% -------------------------------------------------------------------------
function Bezier = bezier(X)
    
    N_VARS      = numel(X);         % degree of freedom
    N_POINTS    = (N_VARS+2)/2;     % number of hinged points for the bezier curve
    DOMAIN      = 0:0.002:1;        % domain for the bezier function
    
    %% generate points for bezier curve    
    points = [0.2 X(1)];            % first point
    
    for i = 2:2:N_VARS-2
        points = cat(1, points, [X(i) X(i+1)]);
    end    
    
    points = cat(1, points, [1 X(N_VARS)]); % last point
   
    %% determining range of bezier over selected domain
    Sigma           = zeros(N_POINTS,1);
    PointMatrix     = zeros(N_POINTS,1);
    BezierMatrix    = [];
    
    % Sigma = (x!/(y!(x-y)!))
    for i=1:N_POINTS
        Sigma(i)=factorial(N_POINTS-1)/(factorial(i-1)*factorial(N_POINTS-i));   
    end

    for i=DOMAIN
        for j=1:N_POINTS
            PointMatrix(j)=Sigma(j)*((1-i)^(N_POINTS-j))*(i^(j-1));
        end

        BezierMatrix=cat(1,BezierMatrix,PointMatrix'); 
    end

    Bezier = BezierMatrix*points;
    
%     line(Bezier(:,1),Bezier(:,2), 'DisplayName', 'Bezier');
%     hold on;
%     line(Points(:,1),Points(:,2), 'DisplayName', 'Line')
%     legend('show');

end
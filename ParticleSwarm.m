% -------------------------------------------------------------------------
% The objective of this function is to :
% determine optimal blade design parameters using Differential Evolution
%
% REFERENCES
% --------------
% Ref1: https://www.mathworks.com/examples/global-optimization/community/20166-particle-swarm-optimization-simulation
%
% INPUTS
% --------------
% lambda        design tip speed ratio of the rotor     [-]
% pitch         design default blade pitch              [deg]
% XC_INIT       inital chord length coefficients        [aN ... a2 a1 a0]
% XT_INIT       initial twist angle coefficients        [aN ... a2 a1 a0]
% profile       profile for chord & twist distribution  [string]
% isPrandtl     is Prandtl correction active            [0/1]
% isGlauert     is Glauert correction active            [0/1]
%
% OUTPUT
% ---------------
% BestDesign    optimal design parameters               [ParticleClass]
%
% CHANGE LOGS
% ---------------
% 08 Sep 2017   initial design
% -------------------------------------------------------------------------

function BestDesign = ParticleSwarm(lambda, pitch, XC_INIT, XT_INIT, profile, isPrandtl, isGlauert)

    disp(sprintf('TSR   \t Loop  \t Cp,max')); %#ok<*DSPS>
    disp(sprintf('----- \t ----- \t -----'));
    t_ps = tic;
    
    %% optimizer parameters
    PS_SWARM_SIZE   = 20;   % number of particles in Particle-Swarm
    PS_INERTIA      = 1.5;  % inertia of the particles in Particle-Swarm
    PS_SELF_ADJUST  = 2.0;  % factor to speed up or slow down the particle based on particle's performance
    PS_SOCIAL_ADJST = 2.0;  % factor to correct the course of particle towards the best particle
    PS_TIME         = 0.2;  % time over which the particles in Particle-Swarm are considered
    PS_N_ITR        = 100;  % number of iterations for Particle-Swarm optimisation    

    %% initialisation
    Swarm(PS_SWARM_SIZE, 1) = ParticleClass();
    for i_particle = 1:PS_SWARM_SIZE
        Particle            = ParticleClass(i_particle/PS_SWARM_SIZE, XC_INIT, XT_INIT);
        Rotor               = BEMRotor(lambda, pitch, Particle.Xc, Particle.Xt, profile, isPrandtl, isGlauert);
        Particle.cP         = Rotor.cP;
        Swarm(i_particle)   = Particle;
    end
    
    %% iteration
    for i_itr = 1:PS_N_ITR            

        %% find out the best particle
        [best_cP, best_index] = max([Swarm(:).cP]);
        BestParticle = Swarm(best_index);

        %% update velocity vectors
        for i_particle = 1:PS_SWARM_SIZE
            Particle = Swarm(i_particle);

            % chord velocity
            c1 = PS_INERTIA*Particle.Xc_velocity;
            c2 = PS_SELF_ADJUST*(Particle.Xc_best - Particle.Xc);
            c3 = PS_SOCIAL_ADJST*(BestParticle.Xc_best - Particle.Xc);

            Particle.Xc_velocity = ((diag(rand(numel(XC_INIT),1))*c1')' + ...
                                   (diag(rand(numel(XC_INIT),1))*c2')' + ...
                                   (diag(rand(numel(XC_INIT),1))*c3')');

            % twist velocity
            t1 = PS_INERTIA*Particle.Xt_velocity;
            t2 = PS_SELF_ADJUST*(Particle.Xt_best - Particle.Xt);
            t3 = PS_SOCIAL_ADJST*(BestParticle.Xt_best - Particle.Xt);

            Particle.Xt_velocity = ((diag(rand(numel(XT_INIT),1))*t1')' + ...
                                      (diag(rand(numel(XT_INIT),1))*t2')' + ...
                                      (diag(rand(numel(XT_INIT),1))*t3')');                

            Swarm(i_particle) = Particle;                    
        end  % end of velocity updation loop    
        
        %% evaluate the position and quality
        for i_particle = 1:PS_SWARM_SIZE 

            Particle    = Swarm(i_particle);
            Particle.Xc = Particle.Xc + PS_TIME*Particle.Xc_velocity;
            Particle.Xt = Particle.Xt + PS_TIME*Particle.Xt_velocity;        

            % check the value of objective function 
            Rotor   = BEMRotor(lambda, pitch, Particle.Xc, Particle.Xt, profile, isPrandtl, isGlauert); 

            if(Rotor.cP > Particle.cP)
                Particle.Xc_best = Particle.Xc;
                Particle.Xt_best = Particle.Xt;
                Particle.cP = Rotor.cP;
            end               

            Swarm(i_particle) = Particle;

        end % end of evaluation

        if(mod(i_itr,10) == 0)
            disp(sprintf('%.1f \t %d \t %.4f', lambda, i_itr, best_cP));
        end
    end % end of iteration  
    
    
    %% return the best design
    t_elapsed       = toc(t_ps);
    BestParticle.t  = t_elapsed;
    BestDesign      = BestParticle;
    disp(sprintf('PSO completed in %.1f mins.\n', t_elapsed/60));

end

% -------------------------------------------------------------------------
% The objective of this function is to :
% determine optimal blade design parameters using Differential Evolution
%
% REFERENCES
% --------------
% Ref1: https://en.wikipedia.org/wiki/Differential_evolution
% Ref2: https://ocw.mit.edu/courses/sloan-school-of-management/15-099-readings-in-optimization-fall-2003/lecture-notes/ses2_storn_price.pdf
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
% BestDesign    optimal design parameters               [IndividualClass]
%
% CHANGE LOGS
% ---------------
% 12 Sep 2017   initial design
% -------------------------------------------------------------------------

function BestDesign = DifferentialEvolution(lambda, pitch, XC_INIT, XT_INIT, profile, isPrandtl, isGlauert)

    disp(sprintf('TSR   \t Loop  \t Cp,max')); %#ok<*DSPS>
    disp(sprintf('----- \t ----- \t -----'));
    t_de = tic;

    %% optimizer parameters
    DE_POPULATION   = 20;                  % number of individuals in the population
    DE_DIFF_WEIGHT  = 0.5;                 % differential weight    
    DE_N_ITR        = 100;                 % number of iterations
    DE_CO_PROB      = 1;                   % cross-over probability
    
    %% initialisation
    Population(DE_POPULATION, 1) = IndividualClass();
    for i_pop = 1:DE_POPULATION
        Individual          = IndividualClass(i_pop/DE_POPULATION, XC_INIT, XT_INIT);
        Rotor               = BEMRotor(lambda, pitch, Individual.Xc, Individual.Xt, profile, isPrandtl, isGlauert);
        Individual.cP       = Rotor.cP;
        Population(i_pop)   = Individual;
    end
    
    %% iteration
    for i_itr = 1:DE_N_ITR
        
        NewPopulation = Population; % next generation
        
        for i_pop = 1:DE_POPULATION
            
            Individual   = Population(i_pop);
            COIndividual = Individual;  % cross-over individual
            
            %% determine agents
            [~, i_best] = max([Population(:).cP]);
            i_agent = [1 1 i_pop];
            
            % each agent and the individual should have unique indices
            while (numel(unique(i_agent)) ~= 3)
               i_agent(1) = randi(DE_POPULATION); 
               i_agent(2) = randi(DE_POPULATION); 
            end
            
            Agent1 = Population(i_agent(1));
            Agent2 = Population(i_agent(2));
            AgentBest = Population(i_best);
            
            %% perform cross-over
            % allow crossover with probability = DE_CO_PROB
            % but at least 1 crossover per individual
            % target-to-best mutation
            
            % chord
            co_point = randi(numel(XC_INIT));            
            for j = 1:numel(XC_INIT)
                if (rand() < DE_CO_PROB || j == co_point)
                    COIndividual.Xc(j) = Individual.Xc(j) + DE_DIFF_WEIGHT*(AgentBest.Xc(j) - Individual.Xc(j)) + ...
                                            DE_DIFF_WEIGHT*(Agent1.Xc(j) - Agent2.Xc(j));
                end
            end
            
            % twist
            co_point = randi(numel(XT_INIT));            
            for j = 1:numel(XT_INIT)
                if (rand() < DE_CO_PROB || j == co_point)
                    COIndividual.Xt(j) = Individual.Xt(j) + DE_DIFF_WEIGHT*(AgentBest.Xt(j) - Individual.Xt(j)) + ...
                                            DE_DIFF_WEIGHT*(Agent1.Xt(j) - Agent2.Xt(j));
                end
            end
            
            
            
            %% check performance of the mutant              
            Rotor               = BEMRotor(lambda, pitch, COIndividual.Xc, COIndividual.Xt, profile, isPrandtl, isGlauert);
            COIndividual.cP     = Rotor.cP;
            
            %% survival of the fittest
            if(COIndividual.cP > Individual.cP)
                NewPopulation(i_pop)   = COIndividual; % the mutant will prevail in the next generation      
            end
                                
        end % end of population
        
        Population = NewPopulation;
        
        %% find out the fittest individual in the population
        [best_cP, best_index] = max([Population(:).cP]);
        
        if(mod(i_itr,10)==0)
            disp(sprintf('%.1f \t %d \t %.4f', lambda, i_itr, best_cP));
        end                        
    end % end of iteration
    
    
    %% return the best design
    t_elapsed       = toc(t_de);    
    BestDesign      = Population(best_index);
    BestDesign.t    = t_elapsed;
    disp(sprintf('DEO completed in %.1f mins.\n', t_elapsed/60));

end % end of function


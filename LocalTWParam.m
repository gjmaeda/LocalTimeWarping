classdef LocalTWParam
    
    properties
        % Initial guess of the time shift. This affects the first parameter of the
        % optimizer theta(1)
        time_shift_guess = 0
        
        % This parameters says that if the cost goes up by more than 
        % stopCriterion_costDeltaUp then the optimization should stop.
        % Cost going up may indicate that the warping of time given by the 
        % optimization is going beyond the limits of the data.       
        stopCriterion_costDeltaUp = 1 %0.2
        
        % Only checks for stopCriterion_costDeltaUp after a certain number
        %  of iterations. This number is specified with costDeltaUpAfterNiter
        stopCriterion_costDeltaUpAfterNiter = 10
        
        % Stop if change of cost is smaller than this number
        stopCriterion_min_rate = 0.005
        
        % limit for optimization (reduced steps is 50)
        iterations  = 1000
        
        % Finite difference perturbation
        perturb = 0.0100
        
        alpha = 0.0100 % Update rate
        nW = 10 % Number of weights on the warping function
        max_number_plots = 0 % Plot figures within the optimization (for debuggin) 
        
        plot_at_iter_number
        h
                
    end
    
    methods
        function obj = LocalDTWParam
            obj.h.scaledCut = 1;
            obj.h.main = 2;
            obj.h.sol = 3;
            obj.h.solpnt = [];
            obj.h.mainpnt = [];
            obj.h.sclpnt = [];
        end
    end
    
end


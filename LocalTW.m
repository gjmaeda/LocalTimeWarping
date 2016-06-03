classdef LocalTW < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        param
        hist
        tref   % reference time
        text   % extended time
        yref
        Phi    % the segment of Phiext that addresses tref
        Phiext % the basis functions on the extended time text
        warping_func
        t_warp
    end
    
    methods
        
        function obj = LocalTW(yref, h)   
            
            obj.param = LocalTWParam();
            
            if isempty(h)
                obj.param.h.scaledCut = [];
                obj.param.h.main = [];
                obj.param.h.sol  = [];                
            else
                obj.param.h.scaledCut = h.scaledCut;
                obj.param.h.main = h.main;
                obj.param.h.sol  = h.sol;
            end
            % Create a fake, normalized time vector for the given reference
            % trajectory
            obj.tref = linspace(0, 1, numel(yref));
            
            % Confirm yref is a column vector.
            [m, n] =  size(yref);
            if m > 1
                obj.yref = yref';
            else
                obj.yref = yref;
            end
            
        end        
        
        function [hist]  = optimize(obj, traj2)

            [m, n] = size(traj2);
            if m > 1
                traj2 = traj2';
            end

            % Resample traj2 to match original size
            % =================================
            if numel(traj2) ~= numel(obj.tref)
                traj2 = interp1(   obj.t01(traj2), traj2, obj.t01(obj.tref) );
            end
            
            % scale data here
            % =================================
            [traj1_, traj2_, ~, r2] = obj.scale(obj.tref, obj.yref, traj2);

            if obj.param.h.scaledCut
                figure(obj.param.h.scaledCut); 
                plot(obj.tref, traj1_, sty('b', [], 4  )  );
                plot(obj.tref, traj2_, sty([0.8 0.8 0.8], [], 4  )  );
                drawnow;
            end
            

            % run main optimization
            % =================================            
            [hist] = obj.gradient_descent(obj.tref, traj1_, traj2_);
            hist.yquery = traj2;
            
            obj.Phi = interp1( obj.text', obj.Phiext', obj.tref')';
            
            % This object contains all data necessary to align other values in
            % relation to the optimized theta
            theta = hist.theta(:,end);
            obj.warping_func = theta(2:end)'*obj.Phiext;        
            obj.t_warp = theta(1) + obj.warping_func.*obj.text;               
            
            %align = Align( hist.theta(:,end), obj.Phiext, obj.text, obj.tref);
            
            hist.pshow = round(linspace(1, hist.k(end),20));
            if hist.pshow(end) ~= hist.k(end)
                hist.pshow(end+1)=hist.k(end);
            end
            
            % recover scale
            % =================================             
            hist = obj.rescale(r2, hist);

            % recover the shift and scale
            if obj.param.h.main
                figure(obj.param.h.main);
                plot(obj.tref, hist.q_sol, sty([0.8 0.8 0.8], [] , 0.1)  )   ;
                plot(obj.tref, hist.q_sol(end,:), sty( 'r', [] , 4)  )   ;
                plot(obj.tref, obj.yref,  sty( 'b', [] , 2)  )   ;
            end
        end

        function [hist] = gradient_descent(obj, t, y, q)
        
            prm = obj.param;
            prm.plot_at_iter_number = round(linspace(1,prm.iterations, prm.max_number_plots));
            if prm.max_number_plots ~=0
                prm.plot_at_iter_number = [prm.plot_at_iter_number  prm.iterations]; 
            end

            % extend time to nExt before and after the movement
            te = obj.get_extended_time(t);
            obj.text = te;
    
            % basis are created for the extended time
            Phi =  obj.create_basis( prm.nW , length(te) );

            % Initial guess of theta
            theta = ones(prm.nW,1);

            % add another theta(1)=0 to represent the shift in time
            theta = [prm.time_shift_guess; theta];

            % extend the values of y and q
            ye = obj.extend_val(t, y, te);
            qe = obj.extend_val(t, q, te);

            % to compute cost function prepare both extended and normal values
            t2.nrm = t;      t2.ext = te;
            q2.nrm = q;      q2.ext = qe;    
            y2.nrm = y;      y2.ext = ye;  

            %% main loop starts here
            hist.k         = [];
            hist.cost      = [];
            hist.q_sol_s   = []; % scaled solution
            hist.theta     = [];
            hist.warpingFunc = [];

            flagStop = 0;     k=1;
            while (k<=prm.iterations) && (flagStop==0)

                [theta_, cost_, q2_, wrpFunc_] = obj.update_theta( t2, y2, q2, Phi, theta, prm, k);
                hist.cost      = [hist.cost cost_];
                hist.theta     = [hist.theta theta_];
                hist.q_sol_s   = [hist.q_sol_s;  q2_];
                hist.k         = [hist.k         k];
                %hist.warpingFunc = [hist.warpingFunc ;  wrpFunc_];
                theta = theta_;
                
                if ~mod(k,10)
                    fprintf('Cost %g\n', cost_ );
                end

               % check if should stop
               if k >= prm.stopCriterion_costDeltaUpAfterNiter
                   if (hist.cost(end)-hist.cost(end-1)) > prm.stopCriterion_costDeltaUp
                       flagStop=1;
                       hist.stop_cause = 'cost went up';
                       disp('  ** Stopping due to cost going up **');
                   end
                   average_rate = mean(abs( diff(hist.cost(end-9:end))) );
                   if prm.stopCriterion_min_rate > average_rate 
                       flagStop=1;
                       hist.stop_cause = 'rate too slow';
                       disp('   ** Stopping due slow convergence **');
                   end
               end
               k=k+1;
            end
            obj.Phiext = Phi;

        end

        function [theta_new, cost, q_aligned, warping_func] = update_theta(obj, t, y, q, Phi, theta, prm, iterNum)

            % compute Jacobian
            % =======================================================================
            [cost, Jacbn, q_aligned, warping_func] = obj.data_jac( t, y, q, Phi, theta, prm, iterNum);

            % Update the parameter
            % ====================================================================
            theta_new = theta - prm.alpha*pinv(Jacbn)*cost';

        end
        
        function [c_0, Jacbn, q_aligned, warping_func] = data_jac(obj, t, y, q, Phi, theta, prm, iterNum)

            % Cost of current parameters
            if sum(prm.plot_at_iter_number == iterNum)
                color_ = 'b';
            else % skip plot
                color_ = [];
            end
            % unperturbed cost
            [c_0, q_aligned, warping_func] = obj.error_cost_function( t, y, q,  theta,  color_, iterNum, Phi ); 

            % compute the partial derivatives
            Jacbn = [];
            for v = 1:length(theta)
                perturb    = zeros(length(theta),1);
                perturb(v) = prm.perturb;
                c_1   = obj.error_cost_function( t, y, q, (theta + perturb) , [], iterNum, Phi); 
                Jp_   = obj.getJacobian(c_1, c_0, prm.perturb);
                Jacbn = [Jacbn  Jp_];
            end

        end
        
        function y_align = unwarp_dofs(obj, y_orig )
            
            [m,n] = size(y_orig);
            if m > n
                y_orig = y_orig';
            end           
            if numel(y_orig(1,:)) ~= numel(obj.tref)
                y_orig = interp1( obj.t01(y_orig(1,:)), y_orig', obj.t01(obj.tref) )';
            end  
            [m,n] = size(y_orig);
            if m > n
                y_orig = y_orig';
            end
            
            y_align =[];
            for d=1:numel(y_orig(:,1))
                % extend the values of y and q
                y_e = obj.extend_val(obj.tref, y_orig(d,:), obj.text);  

                % find the values of the warped q that correspond to the reference time
                y_align = [y_align;  align_data( obj.t_warp, y_e, obj.tref )];
            end
        end
        
        function [] = plot_history(obj, hist)
            
            % get PSI back to the original time
            theta1    = hist.theta(1,:);     % shift in time
            theta2to9 = hist.theta(2:end,:); % weights for bases
            nIter     = numel(hist.theta(1,:));
            
            figurew
            plot(obj.text, obj.Phiext', sty('k', [], 3)  )
            plot(obj.tref, obj.Phi', sty('r', [], 3)  )

            figurew('final_result');
            fWarp2 = [ obj.Phi'*theta2to9 ]'; % warping function
            plot(obj.tref, hist.yquery, sty( [0.6 0.6 0.6], [], 4));
            for k=1:nIter % add the time shift
                t_ref_new(k,:)  = theta1(1,k) + fWarp2(k,:).*obj.tref;
                q_aligned(k,:)  = interp1(  t_ref_new(k,:), hist.yquery, obj.tref );
                if ~mod(k,10) % plot after each N instances
                    plot(obj.tref, q_aligned(k,:), sty_nl([0.8 0.8 0.8], [] , 1));
                end
            end
            plot(obj.tref, obj.yref,  sty('b', [] , 4));
            plot(obj.tref, q_aligned(end,:), sty('r', [] , 4));
            legend({'To\_be\_aligned', 'Reference', 'Aligned\_trajectory'});
            
        end          
        
    end % methods
    
    methods (Static)
 
        function [cost, q_aligned, warping_func] = error_cost_function( t, y, q, theta, style, indx, Phi )

            warping_func = theta(2:end)'*Phi;

            % new time but still extended, you should not compute the cost based on
            % this extended time, but this will avoid interpolating on non existing
            % data
            t.ext_new = theta(1) + warping_func.*t.ext;

            ind = diff(t.ext_new) < 0;
            if sum(ind)
                % What are the consequences of this condition???
                disp('Time is negative');
            end

            % interpolate to find the new values of ob.y in relation to the 
            % timing of the reference (ref.t).
            % find the values of the warped q that correspond to the reference time
            q_aligned = align_data( t.ext_new, q.ext, t.nrm );        

            % Check this condition. It should never happen
            invalid = isnan(q_aligned);
            if sum(invalid)
                error('Interpolation failed');
             end

            N = numel(t.nrm);

            % cost can be computed by direct comparison of y_q2 - y thanks to
            % interpolation
            cost = sum( sqrt( (q_aligned - y.nrm).^2)  )/N * 100 ;
            %fprintf('Cost %g\n', cost );
            
            if ~isempty(style)
                h = figurew(['debug_' num2str(indx)]);
                plot(t.nrm, y.nrm, SGRAYBALL(10));        
                plot(t.nrm, q_aligned, struct('Color', style, 'Marker', 'o'));
                title(['Cost: ' num2str(cost) ' nSize:' num2str(N)]);
            end

        end


        function [J] = getJacobian(c1, c0, perturb)
            dx = c1-c0;
            % comput the Jacbn from the perturbed simulation
            dg1_dp = dx / perturb;
            J = [dg1_dp]';
        end
    
        
        
        function Phi = create_basis(nBasis, nTraj)
            x = linspace(0,1,nTraj);
            sigma= 0.10;
            C=linspace(0,x(end),nBasis);

            Phi=exp(-0.5*(bsxfun(@minus,repmat(x,nBasis,1),C').^2/sigma^2));
            Phi=bsxfun(@rdivide,Phi,sum(Phi));        
            % figurew 'basis'
            % plot(x, Phi, 'b');             
        end
        
        function  ye = extend_val(t, y, te)
            iStart = find( te == t(1) );
            iFinal = find( te == t(end) );
            
            yStart = y(1).*ones(1,   iStart-1);
            yEnd   = y(end).*ones(1, numel(te)-(iFinal)  );            
            ye = [yStart y yEnd];            
        end        
        
        function  te = get_extended_time(t)
            dt = t(2)-t(1);
            nExt = round(numel(t)*0.75);
            tendPlus   = t(end)+dt:dt:(numel(t)+nExt)*dt;
            tstartPlus = -nExt*dt:dt:t(1)-dt;
            te = [ tstartPlus t tendPlus ];
        end
        
        function [traj1_, traj2_, r1, r2] = scale( t, traj1, traj2)

            % to scale try to ignore the beginning and end of demonstrations as
            % they tend to have lots of noise
            i_margin = 1;
            i_r      = i_margin:length(t)-i_margin ;

            % normalize
            r1.nrm.max = max(traj1(i_r));
            r1.nrm.min = min(traj1(i_r));

            r2.nrm.max = max(traj2(i_r));        
            r2.nrm.min = min(traj2(i_r));

            traj1 = traj1./( r1.nrm.max-r1.nrm.min );
            traj2 = traj2./( r2.nrm.max-r2.nrm.min  );

            % shift such that the mean of the trajectory is zero
            max1 = max(traj1(i_r));
            min1 = min(traj1(i_r));
            r1.shift.val = 0.5*( max1+min1 );
            traj1_ = traj1-r1.shift.val;


            max2 = max(traj2(i_r));        
            min2 = min(traj2(i_r));
            r2.shift.val = 0.5*( max2+min2 );
            traj2_ = traj2-r2.shift.val;
        end

        function [hist] = rescale( r, hist)

            rescaled = [];
            nIter = length(hist.pshow);

            for k=1:nIter
                ys  = hist.q_sol_s(hist.pshow(k), :); %scaled
                yr =  (ys + r.shift.val).*(r.nrm.max-r.nrm.min);  %rescaled
                rescaled = [rescaled; yr];
            end

            hist.q_sol = rescaled;

        end
        
        function t = t01(q)
        % return a vector linspace(0,1,N) where N is the size of q
        % 
        % INPUT
        %   q: should be a vector of length N
            [m,n] = size(q);
            if m >1 && n > 1
                error('Vector must be either [Nx1] or [1xN]');
            end
            t = linspace(0,1,numel(q));
            if m > 1
                t = t';
            end
        end       
    end
    
end

function q_aligned = align_data( t_warped, q, t_ref )
% Given 
%    t_warped  as the new warped time
%    t_ref     as the reference time
%    q         as the query trajectroy
%
% return y_aligned
%    as the new values of q that correspond to the reference time

    q_aligned = interp1( t_warped, q, t_ref );
    
end











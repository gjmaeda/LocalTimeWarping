% This code runs a Local Time Warping optimization. It enforces that solutions are
% smooth and do not skip the time indexes as in conventional Dinamic Time Warping [Sakoe, Chiba, 1978].
% Here I assume that enforcing a 1:1 correspondence between indexes is a 
% more natural fit to time warping on dynamical systems since what actually
% differs between two trajectories executed by the robot is the interval
% between the time in indexes. But not the indexes itself. DTW my mess up
% with the indexes, creating unnaturally discontinuous solutions.
%
% LTW may be also useful for optimizing incomplete trajectories,
% although the example here provided does not address it. To my knowledge
% there are very few algorithms that can do DTW on partial trajectories 
% (Dixon, IJCAI, 2005). LTW may be better suited for this task but
% more work is needed.
%
clear; close all; clear; clc; dbstop if error;
addpath('./func_aux'); 

% ==============================================
% Create some data set
% =============================================
if 0
    figurew;
    t = 0:0.05:(2*pi);
    y = 2 + 5.*sin( t ) + 2*sin( t/5 );
    data{1} = y';
    for k=1:5
        y2 = 2+2*abs(randn)*cos(t).*y;
        ts = cumsum(abs( sin(t*5*rand)  + cos(t*5*rand)    ));
        y2 = interp1(t./t(end), y2, ts./ts(end));
        y2 = [ y2(1).*ones(1, round(numel(y)*abs(0.4*randn))  )   y2   y2(end).*ones(1, round(numel(y2)*abs(0.4*randn)) )];    
        plot(y2);
        data{end+1} = y2';
    end
    y = [ y(1).*ones(1, round(numel(y)*0.25)  )   y   y(end).*ones(1, round(numel(y)*0.5) )];
    plot(y, 'r');
    data{1} = y';
    save('data_c.mat', 'data')
else
    load('data_c.mat');
end



nMaxDemos = numel(data);
% ==============================================
% plot misaligned data
% =============================================
hd = figurew('misaligned_data'); 
set_fig_position([0.561 0.0972 0.195 0.792]);
 ylabel 'x'; xlabel 'time steps';
for k=1:nMaxDemos
   plot( linspace(0,1, numel(data{k})), data{k}, sty([0.7 0.7 0.7])  );
end
drawnow;
% Note that here all demonstrations have the same number of time steps.
% This is not a requirement, but the method will resample all solutions to
% the same number of time steps of the reference trajectory.

% ==============================================
% select which training data is the reference and plot it
% =============================================
referenceDemo = 1;  % demonstration used as a reference. 
plot( linspace(0,1,numel(data{referenceDemo}))  , data{referenceDemo}, sty([0.5 0.5 1], [], 3)  );
drawnow;
% Chose which DoF to use as reference [1:x, 2:y, 3:z] for alignment. The other
% DoFs will be aligned according to the solution found for the reference
% DoF. For better results, I recommend you to plot all data and visually 
% inspect which DoF has a consistent, observable, repetitive behavior.
referenceAxis = 1;
refTrajectory = data{referenceDemo}(:, referenceAxis);
refTrajectory = interp1(linspace(0,1,numel(refTrajectory)), refTrajectory, linspace(0,1, 200));

% ==============================================
% Align each of the demonstrations
% =============================================
ltw = LocalTW(refTrajectory', []); % create the object
ltw.param.stopCriterion_costDeltaUp = 1;
for k=1:nMaxDemos
    fprintf('Aligning demonstration %g of %g\n', k, nMaxDemos );
    if k~= referenceDemo
        queryTraj = data{k}(:, referenceAxis);        
        ltw.optimize(queryTraj'); % run the main optimization here
        datan{k} = ltw.unwarp_dofs(  data{k}'  )'; % align all DoFs        
       % plot_result(data{referenceDemo}, data{k}, datan{k}, k );
    else
        datan{k} = data{k};
    end
end

% ==============================================
% Plot the final result
% =============================================
hd = figurew('Aligned_data'); 
set_fig_position([0.761 0.0991 0.195 0.792]);
ylabel 'x'; xlabel 'time steps';
for k=1:nMaxDemos
    plot( linspace(0,1, numel(datan{k})) , datan{k}, sty([0.7 0.7 0.7])  );
end
plot(linspace(0,1, numel(data{referenceDemo})) , data{referenceDemo}, sty([0.5 0.5 1], [], 3)  );










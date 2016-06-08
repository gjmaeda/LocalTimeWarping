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
% The data is real human point-to-point (that is, simple case) with 3 DoFs. 
% Only one DoF is selected to align, and the others are aligned based 
% on the found solution. 
%
clear; close all; clear; clc; dbstop if error;
addpath('./func_aux'); 

% ==============================================
% Load a data set example containing real human demonstrations
% =============================================
load('data2.mat');
nMaxDemos = 7; % limit the number of data to be aligned
randVec = randsample( numel(data), nMaxDemos); % pick some demonstrations at random
data_ = [];
for k = randVec'
    data_{end+1} = data{k};
end
data = data_;

% ==============================================
% plot misaligned data
% =============================================
hd = figurew('misaligned_data'); 
set_fig_position([0.561 0.0972 0.195 0.792]);
subplot(3,1,1); grid on; hold on; ylabel 'x';
subplot(3,1,2); grid on; hold on; ylabel 'y';
subplot(3,1,3); grid on; hold on; ylabel 'z'; xlabel 'time steps';
for k=1:nMaxDemos
    for j=1:3
        subplot(3,1,j);
        plot( linspace(0,1, numel(data{k}(:,j))), data{k}(:,j), sty([0.7 0.7 0.7])  );
    end
end
drawnow;
% Note that here all demonstrations have the same number of time steps.
% This is not a requirement, but the method will resample all solutions to
% the same number of time steps of the reference trajectory.

% ==============================================
% select which training data is the reference and plot it
% =============================================
referenceDemo = 3;  % demonstration used as a reference. 
for j=1:3
    subplot(3,1,j);
    plot( linspace(0,1,numel(data{referenceDemo}(:,j)))  , data{referenceDemo}(:,j), sty([0.5 0.5 1], [], 3)  );
end
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
for k=1:nMaxDemos
    fprintf('Aligning demonstration %g of %g\n', k, nMaxDemos );
    if k~= referenceDemo
        queryTraj = data{k}(:, referenceAxis);        
        ltw.optimize(queryTraj'); % run the main optimization here
        datan{k} = ltw.unwarp_dofs( data{k}'  )'; % align all DoFs        
        plot_result(data{referenceDemo}, data{k}, datan{k}, k );
    else
        datan{k} = data{k};
    end
end

% ==============================================
% Plot the final result
% =============================================
hd = figurew('Aligned_data'); 
set_fig_position([0.761 0.0991 0.195 0.792]);
subplot(3,1,1); grid on; hold on; ylabel 'x';
subplot(3,1,2); grid on; hold on; ylabel 'y';
subplot(3,1,3); grid on; hold on; ylabel 'z'; xlabel 'time steps';
for k=1:nMaxDemos
    for j=1:3
        subplot(3,1,j);
        plot( linspace(0,1, numel(datan{k}(:,j))) , datan{k}(:,j), sty([0.7 0.7 0.7])  );
    end
end
for j=1:3
    subplot(3,1,j);
    plot(linspace(0,1, numel(data{k}(:,j))) , data{referenceDemo}(:,j), sty([0.5 0.5 1], [], 3)  );
end










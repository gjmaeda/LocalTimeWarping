% This is a test code.
%
%
%
%
% This code runs a local DTW optimization. It enforces that solutions are
% smooth and do not skip the time indexes as in conventional DTW.
% It may be also useful for optimizing over incomplete trajectories,
% although the example here provided does not address this part.
%
% Note that although the different demonstrations are not aligned in time,
% the method assumes that all demonstrations are interpolated or resampled 
% such that they all have the same number of time steps.
clear; close all; clear; clc; dbstop if error;
addpath('./func_aux');

% ==============================================
% Load a data set example containing human demonstrations
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
        plot(data{k}(:,j), sty([0.7 0.7 0.7])  );
    end
end
drawnow;

% ==============================================
% select which training data is the reference and plot it
% =============================================
referenceDemo = 3;  % demonstration used as a reference. 
for j=1:3
    subplot(3,1,j);
    plot(data{referenceDemo}(:,j), sty([0.5 0.5 1], [], 3)  );
end
% Chose which DoF to use as reference [1:x, 2:y, 3:z] for alignment. The other
% DoFs will be aligned according to the solution found for the reference
% DoF. It is better to plot all data and visually inspect which DoF has a
% consistent, observable, repetitive behavior. 
referenceAxis = 1;
refTrajectory = data{referenceDemo}(:, referenceAxis);

% ==============================================
% Align each of the demonstrations
% =============================================
dtw = LocalDTW(refTrajectory', []); % create the object
for k=1:nMaxDemos
    fprintf('Aligning demonstration %g of %g\n', k, nMaxDemos );
    if k~= referenceDemo
        queryTraj = data{k}(:, referenceAxis);
        
        dtw.optimize(queryTraj'); % run the main optimization here
        datan{k} = dtw.unwarp_dofs( data{k}'  )'; % align all DoFs
        
        plot_result( data{referenceDemo}, data{k}, datan{k}, k );
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
        plot(datan{k}(:,j), sty([0.7 0.7 0.7])  );
    end
end
for j=1:3
    subplot(3,1,j);
    plot(data{referenceDemo}(:,j), sty([0.5 0.5 1], [], 3)  );
end










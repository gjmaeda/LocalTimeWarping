


clear 
close all
clear 
clc
dbstop if error

addpath('../aux_func/');
addpath('./func_time_warp/');
addpath('./func_new');

% ==============================================
% Load a data set example
% =============================================
load('./2015_02_17_132542_training_raw.mat');

% ==============================================
% put training data in appropriate format for time alignment
% =============================================
exp = Data(training.human, 500);

for k = 1:numel(exp.demos.x(:,1))
    pStart = (0.75).*rand;
    pEnd   = (0.75).*rand;
    x  = fakeQuery( exp.demos.x(k,:), [pStart pEnd] );
    y  = fakeQuery( exp.demos.y(k,:), [pStart pEnd] );
    z  = fakeQuery( exp.demos.z(k,:), [pStart pEnd] );
    
    shiftX = (max(x)-min(x))*(0.35).*rand;
    x = x + shiftX;
    data{k} = [x' y' z'];
end
save('data3.mat', 'data');


break

% ==============================================
% select which training data is the reference and which is to be aligned
% =============================================
referenceTraj = 4;  % demonstration used as a reference. 
demoNumber    = 2;  % demonstration number that you want to align in relation to 
                    % referenceTraj


newRef.x = fakeQuery( exp.demos.x(referenceTraj,:), [0.15 0.015] );
newRef.y = fakeQuery( exp.demos.y(referenceTraj,:), [0.15 0.015] );
newRef.z = fakeQuery( exp.demos.z(referenceTraj,:), [0.15 0.015] );

newQuery.x = fakeQuery( exp.demos.x(demoNumber,:), [0.735 0.255] );
newQuery.y = fakeQuery( exp.demos.y(demoNumber,:), [0.7535 0.255] );
newQuery.z = fakeQuery( exp.demos.z(demoNumber,:), [0.735 0.255] );
                 


prmtmp = start_param();
dtw = LocalDTW(exp.t, newRef.x, prmtmp.h);
[hist] = dtw.optimize(   newQuery.x );
x2 = dtw.unwarp_dofs(newQuery.x);
figurew
plot(newQuery.x, 'bo-')
plot(x2, 'r')

break

dtw.plot_history(hist);




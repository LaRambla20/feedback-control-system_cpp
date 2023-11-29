clear all
close all
clc

% Once the C++ code has been run and the scope_data_mod.tsv file has been saved: 

%% Plot of the reference data
table = readtable("./../data/scope_data.tsv", "FileType","text",'Delimiter', '\t');

time = table2array(table(:,1));
reference = table2array(table(:,2));
plant_reference = table2array(table(:,3));
plant_output = table2array(table(:,4));

figure
plot(time, reference, 'LineWidth', 3.0);
hold on
title('Simulink Signals')
grid on
plot(time, plant_reference, 'LineWidth', 3.0);
plot(time, plant_output, 'LineWidth', 3.0);
legend('inRef','plantRef','plantOut')

%% Plot of the generated data
table = readtable("./../data/scope_data_mod.tsv", "FileType","text",'Delimiter', '\t');

time = table2array(table(:,1));
reference = table2array(table(:,2));
plant_reference = table2array(table(:,3));
plant_output = table2array(table(:,4));

figure
plot(time, reference, 'LineWidth', 3.0);
hold on
title('C++ Signals')
grid on
plot(time, plant_reference, 'LineWidth', 3.0);
plot(time, plant_output, 'LineWidth', 3.0);
legend('inRef','plantRefCpp','plantOutCpp')
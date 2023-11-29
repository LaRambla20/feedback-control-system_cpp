% Once the simulink model has been run: 

%% Save and plot the data output by the simulink simulation
table = array2table([out.simout.time out.simout.data]);
writetable(table,'./../data/scope_data.tsv', 'filetype','text', 'delimiter','\t','WriteVariableNames',false);

time = out.simout.time;
reference = out.simout.data(:,1);
plant_reference = out.simout.data(:,2);
plant_output = out.simout.data(:,3);

figure
plot(time, reference, 'LineWidth', 3.0);
hold on
title('Simulink Signals')
grid on
plot(time, plant_reference, 'LineWidth', 3.0);
plot(time, plant_output, 'LineWidth', 3.0);
legend('inRef','plantRef','plantOut')
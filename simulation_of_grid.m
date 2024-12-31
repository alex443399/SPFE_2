%% Loading the data
Data = readtable('household_power_consumption.txt');
% Obtained from https://archive.ics.uci.edu/dataset/235/individual+household+electric+power+consumption

%%
Power_readings_50_hours = Data(1:3000,7:9);
size(Power_readings_50_hours)
input_signal_power_in_watthours_per_minute = table2array(Power_readings_50_hours)';
% Units according to the link
input_signal_power_in_watts = input_signal_power_in_watthours_per_minute / 60;
%% Model data
% Network components
number_of_generators = 3;
number_of_loads = 3;
number_of_components = number_of_loads + number_of_generators;

% Values in table 1
array_k_governor_feedback_gain = [100 125 150];
array_R_droop_characteristic = [0.05 0.04 0.033];
array_D_rotor_damping = ones(1,6) * 1.5;
array_M_rotor_inertia = [10 5 8.33 1 3.33 3.33];
array_T_machine_time_constant = [5 3 4];

% Suceptance matrix
B = [
    -0.334 0 0 0.176 0.158 0;
    0 -0.455 0 0.306 0 0.149;
    0 0 -0.567 0 0.358 0.209;
    0.176 0.306 0 -0.482 0 0;
    0.158 0 0.358 0 -0.518 0;
    0 0.149 0.209 0 0 -0.358
];

% I think the diagonal is the sum of the rest of the row,
% and the values outisde the diagonal are in Table 1

%% We create the system
Full_order_model_discrete = Get_discrete_time_model(number_of_generators, ...
    number_of_loads, array_k_governor_feedback_gain, ...
    array_R_droop_characteristic, array_D_rotor_damping, ...
    array_M_rotor_inertia, array_T_machine_time_constant, ...
    B);


%% Simulating
y = lsim(Full_order_model_discrete, input_signal_power_in_watts);
%% Plotting results
plot(y())
%% Plot to compare with input
figure
subplot(2,1,1)
plot(y())
subplot(2,1,2)
plot(input_signal_power_in_watts')
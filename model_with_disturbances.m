%% We take the simulation results and obtain the disturbances
% Since P=B\delta, we have that P_{1,5} = B_{1,5} (\delta(5)-\delta(1)), idem for
% P_{2,6} = B_{2,6} (\delta(6)-\delta(2)). This is because the power in an
% interconnection depends on the difference between shaft angles. However,
% by convention \delta(1) = 0. In actuality, our variables \delta_i are
% actually \delta_i-\delta_1. This is because their derivatives are
% \omega_i-\omega_1.

matrix_to_get_disturbances = [
    -0 0 0 0.158 0;
    -0.149 0 0 0 0.149
];

disturbance_power = matrix_to_get_disturbances * x_simulation_state(:,13:end)';

plot(t_span/3600, disturbance_power(1,:), t_noisy/3600, disturbance_power(2,:))
legend('Power from bus 5\rightarrow 1', 'Power from bus 6\rightarrow 2')
title('Power connecting both parts of the network')
ylabel('Power [Watts]')
xlabel('Time [hours]')

%% Model
Partial_order_model_discrete = Get_discrete_time_partial_model(number_of_generators, ...
    number_of_loads, array_k_governor_feedback_gain, ...
    array_R_droop_characteristic, array_D_rotor_damping, ...
    array_M_rotor_inertia, array_T_machine_time_constant, ...
    B)
%% Simulation
[y_partial, ~, x_partial]= lsim(Partial_order_model_discrete, disturbance_power);
plot(y_partial)
%% Validation
plot(y_partial(:,1))
hold on
plot(x_simulation_state(:,3));
hold off
Q = diag([1e-4, 1e-4, 1e-7, 1e-4, 1e-4, 1e-7, 1e-7, 1e-3, 1e-3]);
R = diag([1e-8, 1e-8]);


w = mvnrnd(zeros(9,1), Q, length(x_simulation_state));
v = mvnrnd(zeros(2,1), R, length(x_simulation_state));

u = [disturbance_power;
    w';
    v'];
%%
model_with_noise_channels = Partial_model_with_noise(Partial_order_model_discrete);

[y_noisy, ~, x_noisy] = lsim(model_with_noise_channels, u);

%% Running SISE filter

[d_hat, x_hat] = SISE_filter(y_noisy, Q, ...
    Partial_order_model_discrete.A, ...
    Partial_order_model_discrete.B, ...
    Partial_order_model_discrete.C, Q, R);
%%
y_hat = Partial_order_model_discrete.C * x_hat';
plot(y_hat'-y_noisy)
%%
plot(t_span_hr,x_noisy(:,4),t_span_hr,x_hat(:,4))
legend('state', 'SISE')
xlabel('Time [hours]')
title('valve position bus 2')
%%
plot(t_span_hr,x_noisy(:,5),t_span_hr,x_hat(:,5))
legend('state', 'SISE')
xlabel('Time [hours]')
title('Mechanical power bus 2')
%%
plot(t_span_hr,x_noisy(:,6),t_span_hr,x_hat(:,6))
legend('state', 'SISE')
xlabel('Time [hours]')
title('Angular velocity bus 2')
%%
plot(t_span_hr,x_noisy(:,8),t_span_hr,x_hat(:,8))
legend('state', 'SISE')
xlabel('Time [hours]')
title('Shaft angle bus 2')
%%
plot(t_span_hr,disturbance_power(1,:),t_span_hr,d_hat(:,1))
legend('state', 'SISE')
xlabel('Time [hours]')
title('Power flow from bus 5\rightarrow 1')

%%
function nsys = Partial_model_with_noise(dsys)
    A = dsys.A;
    D = dsys.D;
    B = [dsys.B eye(length(A)) zeros(length(A), length(D))];
    C = dsys.C;
    D = [zeros(size(dsys.D)) zeros(length(D), length(A)) eye(size(D))];
    nsys = ss(A,B,C,D,dsys.Ts);
end


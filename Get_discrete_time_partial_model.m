function sys = Get_discrete_time_partial_model(number_of_generators, ...
    number_of_loads, arr_k, arr_R, arr_D, arr_M, arr_T, mat_B_sucept)
    %% A matrix cont time

    G_1 = [
            -arr_k(1)*arr_R(1) 0 arr_k(1);
            1/arr_T(1) -1/arr_T(1) 0;
            0 -1/arr_M(1) -arr_D(1)/arr_M(1)
            ];
    G_2 = [
            -arr_k(2)*arr_R(2) 0 arr_k(2);
            1/arr_T(2) -1/arr_T(2) 0;
            0 -1/arr_M(2) -arr_D(2)/arr_M(2)
            ];
    L_4 = -arr_D(4)/arr_M(4);
    A_upper_left_corner = blkdiag(G_1,G_2,L_4);
    A_lower_left = [
        zeros(2) -ones(2,1) zeros(2) eye(2)
        ];
    A_upper_right_corner = [
        zeros(2);
        [mat_B_sucept(1,2)/arr_M(1) mat_B_sucept(1,4)/arr_M(1)];
        zeros(2);
        [mat_B_sucept(2,2)/arr_M(2) mat_B_sucept(2,4)/arr_M(2)];
        [mat_B_sucept(4,2)/arr_M(4) mat_B_sucept(4,4)/arr_M(4)];
    ];

    A_cont_time = [
        A_upper_left_corner A_upper_right_corner;
        A_lower_left zeros(2)
        ];
    %% B and C matrices cont time
    B_cont_time = [
        zeros(2);
        [1/arr_M(1) 0];
        zeros(2);
        [0 1/arr_M(2)];
        zeros(3,2);
        ];
    C_cont_time = [
        kron(eye(2), [0 0 1]) zeros(2,3)];
    %% System
    sys_cont = ss(A_cont_time, B_cont_time, C_cont_time, zeros(2));

    % Discretizing
    % For the sake of simplicity we assume zero order hold
    Ts = 60; % Model units are in second, data is measured every minute
    sys = ss(c2d(sys_cont, Ts),'zoh');

end


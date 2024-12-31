function sys = Get_discrete_time_model(number_of_generators, ...
    number_of_loads, arr_k, arr_R, arr_D, arr_M, arr_T, mat_B_sucept)
    % We get the matrix script_A defined in equation (5)
    % To make it we constuct 4 submatrices, it is further explained inside
    % the method
    script_A = get_continous_time_state_matrix( ...
    number_of_generators, number_of_loads, arr_k, arr_R, arr_D, ...
    arr_M, arr_T, mat_B_sucept);

    % We write the input matrix. The input is the load power consumption,
    % the interaction is defined in equation (2). This matrix is consistent
    % with that
    script_B = get_continous_time_input_matrix(arr_M);
    
    % Creating continous system
    % For now we will let C = I_n so we can observe everything, 
    % as we need those to later make the simulations with 
    % limited information.
    n = length(script_A);
    continous_time_system = ss(script_A, script_B, eye(n), 0);

    % Discretizing
    % For the sake of simplicity we assume zero order hold
    Ts = 60; % Model units are in second, data is measured every minute
    sys = ss(c2d(continous_time_system, Ts),'zoh');
    
end


function script_A = get_continous_time_state_matrix( ...
    number_of_generators, number_of_loads, arr_k, arr_R, arr_D, arr_M, arr_T, mat_B_sucept)
    % We create the block matrices for the system matrix
    % First the G and L matrices of (4)
    G = zeros(number_of_generators, 3, 3);
    L = zeros(number_of_loads, 1, 1);
    
    for i = 1:number_of_generators
        k = arr_k(i);
        R = arr_R(i);
        D = arr_D(i);
        M = arr_M(i);
        T = arr_T(i);
        G(i,:,:) = [
            -k*R 0 k;
            1/T -1/T 0;
            0 -1/M -D/M
            ];
    end
    
    for i = 1:number_of_loads
        D = arr_D(i+number_of_generators);
        M = arr_M(i+number_of_generators);
        L(i,:,:) = -D/M;
    end
    
    % We create the top right submatrix
    I_matrix = [
        zeros(2,5);
        mat_B_sucept(1,2:6)/arr_M(1);
        zeros(2,5);
        mat_B_sucept(2,2:6)/arr_M(2);
        zeros(2,5);
        mat_B_sucept(3,2:6)/arr_M(3);
        mat_B_sucept(4,2:6)/arr_M(4);
        mat_B_sucept(5,2:6)/arr_M(5);
        mat_B_sucept(6,2:6)/arr_M(6)
    ];
    
    % Top left submatrix, the block diagonal of the dynamics of the generators and loads 
    submatrix_GL = blkdiag(...
        reshape(G(1,:,:),[3,3]),...
        reshape(G(2,:,:),[3,3]),...
        reshape(G(3,:,:),[3,3]),...
        reshape(L(1),[1,1]),...
        reshape(L(2),[1,1]),...
        reshape(L(3),[1,1]));
    
    % Bottom left submatrix, consisting of 0,+1,-1. We construct the
    % transpose because it is easier
    A_lower_left_submatrix = [
        zeros(2,5);
        -ones(1,5);
        zeros(2,5);
        1 0 0 0 0;
        zeros(2,5);
        0 1 0 0 0;
        0 0 1 0 0;
        0 0 0 1 0;
        0 0 0 0 1;
        ]';
    
    % Matrix script A in (5)
    % This is the state matrix of our system
    script_A = [
        submatrix_GL I_matrix;
        A_lower_left_submatrix zeros(5)];
end

function script_B = get_continous_time_input_matrix(arr_M)
    %% Creating input matrix
    % There are 3 inputs, corresponding to the load powers. These drive the
    % system as stated in Page 3 right before section III starts. These are
    % obtained from the data file.
    % We have in (2) that
    % $\frac{d}{dt} [\omega_r]_i = ... - 1/M_i * P_Li$
    % So that means that the derivatives depend on these powers
    % In x, the indeces of [\omega_r]_i are 3*3+1, 3*3+2, 3*3+3 because they
    % come after the first three components. This means that they are 10,11,12
    
    % Input matrix of our system
    
    script_B = [
        zeros(9,3);
        diag([...
        -1/arr_M(4), ...
        -1/arr_M(5), ...
        -1/arr_M(6)]);
        zeros(5,3)
    ];
end
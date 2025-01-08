function [d_hat, x_hat] = SISE_filter(y_sim, P0, A, G, C, Q, R)
    T_MAX = length(y_sim);
    n = length(A);
    x_hat = zeros(T_MAX, n);
    d_hat = zeros(T_MAX, length(G(1,:)));

    P = P0;
    X = 0;
    K = 0;

    x_hat_current = zeros(n,1);

    for t = 1:T_MAX
        X = A * P * A' + Q;
        S = C * X * C' + R;
        K = X * C' / S;
        M = inv( ...
            (G'*C'/S)*C*G) * G' * C'/S;
        P = (eye(n) - K*C) * ((eye(n) - G*M*C)*X*(eye(n)-G*M*C)' + G * M * R * M' * G');
        y = y_sim(t,:)';
        d_hat(t,:) = M * (y - C * A * x_hat_current);
        
        innovation = K * (y - C * A * x_hat_current - C * G * d_hat(t,:)');

        x_hat_current = A * x_hat_current + G * d_hat(t,:)' + innovation;
        x_hat(t,:) = x_hat_current;
    end
    
end


function [sum_rate, max_r, min_r] = sum_rate_oma(H, Po, P_t)
    % H: M x K channel matrix (M antennas, K users)
    % Po: M x K precoding matrix (each column is a precoder for user k)
    % P_t: transmit power
    
    
    [M, K] = size(H);    % M antennas, K users
    sum_rate_oma = 0;

    R1 = 0;
    R2 = 0;
    
    % Power allocation
    P = sqrt(P_t)*Po;
    alp = 0.5;

    for k = 1:K
        hk = H(:, k);         % Channel vector to user k (M x 1)
        pk = P(:, k);         % Precoding vector for user k (M x 1)
        
        signal_power = abs(hk' * pk)^2;

        Rk = alp * log2(1 + signal_power);

        if k == 1
            R1 = Rk;  % Store the rate for user 1
        else
            R2 = Rk;  % Store the rate for other users
        end

        sum_rate_oma = sum_rate_oma + Rk;
    end
    sum_rate = sum_rate_oma;
    max_r = max(R1, R2);
    min_r = min(R1, R2);
end

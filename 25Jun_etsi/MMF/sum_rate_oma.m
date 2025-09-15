function [sum_rate, max_r, min_r] = sum_rate_oma(H, Po, P_t)
    % H: M x K channel matrix (M antennas, K users)
    % Po: M x K precoding matrix (each column is a precoder for user k)
    % P_t: transmit power
    
    
    [M, K] = size(H);    % M antennas, K users

    R1 = 0;
    R2 = 0;
    R_SUM_MAX = 0;
    
    % Power allocation
    P = sqrt(P_t)*Po;

    ti = [0:0.01:1];

    for i = 1:length(ti)
        
        alp = ti(i);

        sum_rate_oma = 0;

        for k = 1:K
            hk = H(:, k);         % Channel vector to user k (M x 1)
            pk = P(:, k);         % Precoding vector for user k (M x 1)
            
            signal_power = abs(hk' * pk)^2;
    
            Rk = log2(1 + signal_power);
    
            if k == 1
                R1 = (1 - alp) * Rk;  % Store the rate for user 1
            else
                R2 = (alp) * Rk;  % Store the rate for other users
            end
    
            sum_rate_oma = sum_rate_oma + Rk;
        end

        if min(R1, R2) > R_SUM_MAX
            
            R_SUM_MAX = min(R1, R2);

            R_1 = R1;
            R_2 = R2;

        end

    end
    sum_rate = (R_1 + R_2) / 2;
    max_r = max(R_1, R_2);
    min_r = min(R_1, R_2);
end

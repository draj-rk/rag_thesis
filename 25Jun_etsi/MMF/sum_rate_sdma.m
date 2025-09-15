function [sum_rate, max_r, min_r] = sum_rate_sdma(H, Po, P_t)
    % H: M x K channel matrix (M antennas, K users)
    % Po: M x K precoding matrix (each column is a precoder for user k)
    % P_t: transmit power
    
    
    [M, K] = size(H);    % M antennas, K users
    sum_rate_sdma = 0;

    R1 = 0;
    R2 = 0;
    
    % Power allocation
    P = sqrt(P_t/2)*Po;

    for k = 1:K
        hk = H(:, k);         % Channel vector to user k (M x 1)
        pk = P(:, k);         % Precoding vector for user k (M x 1)
        
        signal_power = abs(hk' * pk)^2;

        interference = 0;
        for j = 1:K
            if j ~= k
                pj = P(:, j);
                interference = interference + abs(hk' * pj)^2;
            end
        end

        Rk = log2(1 + (signal_power / (interference + 1)));

        if k == 1
            R1 = Rk;  % Store the rate for user 1
        else
            R2 = Rk;  % Store the rate for other users
        end

        sum_rate_sdma = sum_rate_sdma + Rk;
    end
    sum_rate = sum_rate_sdma / 2;
    max_r = max(R1, R2);
    min_r = min(R1, R2);
end

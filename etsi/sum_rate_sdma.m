function sum_rate_sdma = sum_rate_sdma(H, Po, Pt)
    % H: M x K channel matrix (M antennas, K users)
    % P: M x K precoding matrix (each column is a precoder for user k)
    
    [M, K] = size(H);    % M antennas, K users
    sum_rate_sdma = 0;
    
    % Power allocation
    P = sqrt(Pt/2)*Po;

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

        Rk = log2(1 + signal_power / (interference + 1));
        sum_rate_sdma = sum_rate_sdma + Rk;
    end
end

function [sum_rate, max_r, min_r] = sum_rate_NOMA(H, Po, P_t)
    % H: M x K channel matrix (M antennas, K users)
    % Po: M x K precoding matrix (each column is a precoder for user k)
    % P_t: transmit power
    
    
    [M, K] = size(H);    % M antennas, K users

    R1 = 0;
    R2 = 0;
    R_SUM_MAX = 0;
    
    % Power allocation
    P = sqrt(P_t)*Po;

    % Precoders
    p1 = P(:, 1);
    [U, ~, ~] = svds(H, 1);
    p2 = sqrt(P_t)*U;

    ti = [0:0.01:1];

    for i = 1:length(ti)
        
        tk = ti(i);

        sum_rate_noma = 0;

        % User 1
        R1 = log2(1 + abs(H(:, 1)' * sqrt(tk) * p1)^2);

        for k = 1:K
            hk = H(:, k);         % Channel vector to user k (M x 1)
            %pk = P(:, k);         % Precoding vector for user k (M x 1)
            
            signal_power = abs(hk' * sqrt(1 - tk) * p2)^2;
            inter = 1 + abs(hk' * sqrt(tk)*p1)^2;
    
            Rk = log2(1 + (signal_power/inter));
    
            if k == 1
                R2_1 = Rk;  
            else
                R2_2 = Rk;  
            end
        end

        R2 = min(R2_1, R2_2);
    
        sum_rate_noma = R1 + R2;

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

function rsma_rate = sum_rate_rsma(H, Po, Pt)
    % H: M x K channel matrix (M antennas, K users)
    % Po: M x K precoding matrix (each column is a precoder for user k)
    
    [M, K] = size(H);    % M antennas, K users
    sum_rate_rsma = 0;
    
    R_SUM_MAX = 0;
    
    ti = [0:0.05:1];
    
    for i = 1:length(ti)
        
        t = ti(i);
        
        % Compute the common precoder
        p_c = multicast_common_pre(H, Po);
    
        % Power allocation
        p_c = sqrt(Pt*t) * p_c;
        % Power allocation to private precoders
        P = sqrt(Pt*(1 - t) / 2)*Po;
        

        % For common rates
        Rc_vals = zeros(2, 1);
        for k = 1:K
            hk = H(:, k);
            signal = abs(hk' * p_c)^2;
            inter = 0;

            for j = 1:K
                inter = inter + abs(hk' * P(:,j))^2;
            end

            Rc_vals(k) = log2(1 + signal / (inter + 1));
        end

        Rc = min(Rc_vals);


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
            sum_rate_rsma = sum_rate_rsma + Rk;
        end

        sum_rate_rsma = sum_rate_rsma + Rc;
        
        if sum_rate_rsma > R_SUM_MAX
            
            R_SUM_MAX = sum_rate_rsma;
        end
        
    end
    
    rsma_rate = R_SUM_MAX;
    
end
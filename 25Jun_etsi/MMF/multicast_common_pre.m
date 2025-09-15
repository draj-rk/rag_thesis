function pc = multicast_common_pre(H, Po, P_t, t)
    % H: M x K channel matrix (M antennas, K users)
    % Po: M x K precoding matrix (each column is a precoder for user k)
    % P_t: transmit power
    
    h1 = H(:,1);
    h2 = H(:,2);
    
    beta1_sq = 0;
    beta2_sq = 0;
    
    P = (P_t*(1 - t) / 2)*Po;
    
    for r = 1:2
        beta1_sq = beta1_sq + abs(h1' * P(:,r))^2;
        beta2_sq = beta2_sq + abs(h2' * P(:,r))^2;
    end
    
    beta1_sq = beta1_sq + 1;
    beta2_sq = beta2_sq + 1;
    
    h1_tilde = h1 / sqrt(beta1_sq);
    h2_tilde = h2 / sqrt(beta2_sq);
    
    H_tilde = [h1_tilde h2_tilde];
    
    G = H_tilde' * H_tilde;
    
    a11 = real(G(1,1));
    a22 = real(G(2,2));
    a12 = G(1,2);
    
    abs_a12 = abs(a12);
    denom = a11 + a22 - 2*abs_a12;
    lambda = (a11 * a22 - abs_a12^2) / denom;
    
    mu_vec = (1/denom) * [a22 - abs_a12; a11 - abs_a12];
    
    a12_phase = angle(a12);
    
    pc_bar = (mu_vec(1) * h1_tilde) + (mu_vec(2) * h2_tilde * exp(-1i * a12_phase));
    
    pc = pc_bar / sqrt(lambda);
end
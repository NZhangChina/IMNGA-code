function [H, F_v, F_plus_v, A, A_plus_v, result_imcan1_all, delte_v_all,Tobj] = ...
      Imc_IMNGA(X_miss, Y_miss, A_v, F_v, G_v, alpha1, alpha2, alpha3,  ...
      max_iter_out, num_samples_v,alpha0, k_neigh,isnormalization2)

if nargin < 12
    k_neigh = 100; %{10 15 20 30 35 40 45 50 80 100}
end
if nargin < 13
    isnormalization2 = 0;
end
num_class = max(Y_miss);
num_samples = size(Y_miss,1);
num_view = size(A_v,2);

num_imX = 0; distX = 0; 
for vv = 1:num_view
    P_v{1,vv} = G_v{1,vv}' * ones(size(G_v{1,vv},1),num_class);
    I_v{1,vv} =  eye(num_samples) .* (ones(num_samples,size(G_v{1,vv},1)) * G_v{1,vv});
    temp_dX = L2_distance_1(X_miss{1,vv}', X_miss{1,vv}');
    distX_v{1,vv} = G_v{1,vv}' * temp_dX/mean(mean(temp_dX)) * G_v{1,vv};
    distX = distX + distX_v{1,vv};
    W_v{1,vv} = G_v{1,vv}' * ones(num_samples_v(vv),num_samples_v(vv)) * G_v{1,vv};
    num_imX = num_imX + W_v{1,vv};
end
distX = distX./(num_imX + eps);
[distX1, idx] = sort(distX,2);
A = zeros(num_samples);
rr = zeros(num_samples,1);
% k_neigh = 5;
for i = 1:num_samples
    di = distX1(i,2:k_neigh+2);
    rr(i) = 0.5*(k_neigh*di(k_neigh+1)-sum(di(1:k_neigh)));
    id = idx(i,2:k_neigh+2);
    A(i,id) = (di(k_neigh+1)-di)/(k_neigh*di(k_neigh+1)-sum(di(1:k_neigh))+eps);
end
r = 0; %mean(rr);
A0 = A-diag(diag(A));
A0 = (A0+A0')/2;

for vv = 1:num_view
    temp_A_plus_v = (A + A_v{1,vv})./(W_v{1,vv}+1);
    A_plus_v{1,vv} = zeros(num_samples);
    for i = 1:num_samples
        a0 = temp_A_plus_v(i,:);
        idxa1 = find(a0>0);
        ad = temp_A_plus_v(i,idxa1);
        A_plus_v{1,vv}(i,idxa1) = EProjSimplex_new(ad);
    end

    temp_A0_plus = (A_plus_v{1,vv}+A_plus_v{1,vv}')/2;
    if isnormalization2 == 0
        [F_plus_v{1,vv},~,~] = svds(temp_A0_plus,num_class);
    else
        temp_D10_plus = diag(1./sqrt(sum(temp_A0_plus, 2)+eps));
        temp_A0_plus = temp_D10_plus*temp_A0_plus*temp_D10_plus;
        [F_plus_v{1,vv},~,~] = svds(temp_A0_plus,num_class);
    end
end

iter_inner = 0;
IsConverge_out = 0;
while (IsConverge_out == 0&&iter_inner<max_iter_out)
    iter_inner = iter_inner +1;

    %% Update delte_v 
    if iter_inner == 1
        delte_v = ones(num_view,1)/num_view;
        sum_delte_v = 1;
    else
        sum_delte =0;
        for vv = 1:num_view
            temp1 = A_v{1,vv} - (H.*P_v{1,vv})*F_v{1,vv}';
            temp2 = A_plus_v{1,vv} - H*F_plus_v{1,vv}';
            temp3 = distX_v{1,vv} .* A0;
            temp4 = abs(A_plus_v{1,vv}.*W_v{1,vv} - A_v{1,vv});
            temp5 = abs(A_plus_v{1,vv} - A);
            delte_v(vv,1) = alpha0 * (norm(temp1,'fro'))^2 + (norm(temp2,'fro'))^2 + ... 
                alpha2 * (norm(temp3,'fro'))^2 + alpha3 * (norm(temp4,'fro'))^2 + ...
                alpha3 * (norm(temp5,'fro'))^2;
            sum_delte = sum_delte + delte_v(vv,1);
            Tobj(iter_inner,:) = [(norm(temp1,'fro'))^2 (norm(temp2,'fro'))^2 ...
                (norm(sum_distX_v.*A0,'fro'))^2 (norm(temp3,'fro'))^2 ... 
                (norm(temp4,'fro'))^2 (norm(temp5,'fro'))^2];
        end
        sum_delte_v =0;
        for vv = 1:num_view
            delte_v(vv,1) = (delte_v(vv,1)/sum_delte)^(-1/2);
            sum_delte_v = sum_delte_v + delte_v(vv,1);
        end
    end

    %% Update H 
    I_delte_v = zeros(num_samples,num_samples);
    H = zeros(num_samples,num_class);
    for vv = 1:num_view
        H = H + delte_v(vv,1)*(alpha0 * A_v{1,vv}*F_v{1,vv}+A_plus_v{1,vv}*F_plus_v{1,vv});
        I_delte_v = I_delte_v + delte_v(vv,1) * (alpha0 * I_v{1,vv}+eye(num_samples));
    end
    L_delte_v = diag(sum(A0)) - A0;
    H = (I_delte_v + 2 * alpha1 * L_delte_v + eps) \ H;
    H = max(H, 0);

    %% Update A_v
    for vv = 1:num_view
        A_plus_v{1,vv} = zeros(num_samples,num_samples);
        temp_A_plus_v = A + A_v{1,vv}+ A_plus_v{1,vv};
        dA_plus_v = H*F_plus_v{1,vv}' + alpha3*A + alpha3*A_v{1,vv};
        dA_plus_v = dA_plus_v./(1+alpha3+alpha3*W_v{1,vv});
        for i = 1:num_samples
            a0 = temp_A_plus_v(i,:);
            idxa1 = find(a0>0);
            ad = dA_plus_v(i,idxa1);
            A_plus_v{1,vv}(i,idxa1) = EProjSimplex_new(ad);
        end
    end

    %% Update A 
    distH = L2_distance_1(H',H');
    sum_A_plus_v = 0;
    sum_distX_v = 0;
    for vv = 1:num_view
        sum_A_plus_v = A_plus_v{1,vv}*delte_v(vv,1) + sum_A_plus_v;
        sum_distX_v  =  distX_v{1,vv}*delte_v(vv,1) + sum_distX_v;
    end
    A = zeros(num_samples,num_samples);
    for i=1:num_samples
        idxa1 = idx(i,2:k_neigh+1);
        a0 = sum_A_plus_v(i,:);
        idxa2 = find(a0>0);  % unique(A)
        idxa0 = unique([idxa1 idxa2]);
        dfi = distH(i,idxa0);
        dxi = sum_distX_v(i,idxa0);
        dAvi = sum_A_plus_v(i,idxa0);
        ad = (2*alpha3*dAvi-alpha2*dxi-alpha1*dfi)/(2*r*alpha2+2*sum_delte_v*alpha3);
        A(i,idxa0) = EProjSimplex_new(ad);
    end
    A0 = A-diag(diag(A));
    A0 = (A0+A0')/2;

    %% Update F_v,F_plus_v 可以将这一步放在H后
    for vv = 1:num_view
        temp = A_v{1,vv}'*(H.*P_v{1,vv});
        [Ur, ~, Vr] = svds(temp,num_class);
        F_v{1,vv} = Ur*Vr';
        temp = A_plus_v{1,vv}'*H;
        [Ur, ~, Vr] = svds(temp,num_class);
        F_plus_v{1,vv} = Ur*Vr';
    end    

    %%
    [non ,Y_Hpred] = max(H,[],2);
    % [ACC MIhat Purity]
    temp_res = ClusteringMeasure(Y_miss, Y_Hpred);
    result_imcan1_all(iter_inner,:) = temp_res;
    delte_v_all(iter_inner,:) = delte_v;

%     % convergence checking
%     if iter_inner>1
%         temp_obj = Tobj(iter_inner -1);
%     else
%         temp_obj = 0;
%     end
%     if abs(obj - temp_obj)/temp_obj <1e-8
%         IsConverge_out = 0;
%     end
end

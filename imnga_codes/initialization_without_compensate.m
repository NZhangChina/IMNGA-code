function [A_v, F_v, G_v, M_v, num_samples_v] = initialization_without_compensate(X_miss,...
    X_miss_id,num_class,num_samples,k_neigh_init,issymmetric,isnormalization)


if nargin < 6
    k_neigh_init = 100; %{10 15 20 30 35 40 45 50 80 100}
end
if nargin < 6
    issymmetric = 1;
end
if nargin < 7
    isnormalization = 1;
end
num_view = size(X_miss,2);
num_samples_v = zeros(num_view,1);
for vv = 1:num_view
    temp_X_v_miss = X_miss{1,vv};
    num_samples_v(vv,1) = size(temp_X_v_miss,1);
    
    G_v{1,vv} = zeros(num_samples_v(vv,1),num_samples);
    temp_X_v_id = X_miss_id{1,vv};
    for sample_ii = 1:num_samples_v(vv,1)
        G_v{1,vv}(sample_ii,temp_X_v_id(sample_ii))=1;
    end
    
    if issymmetric == 0 
        M_v{1,vv} = G_v{1,vv}' * ones(num_samples_v(vv,1),num_samples);
    else
        M_v{1,vv} = G_v{1,vv}' * ones(num_samples_v(vv,1),num_samples_v(vv,1)) * G_v{1,vv};
    end
    
    A_v_sub{1,vv} = constructW_PKN(temp_X_v_miss',k_neigh_init, issymmetric);
    A_v_miss{1,vv} = G_v{1,vv}'*A_v_sub{1,vv}*G_v{1,vv};
    
    temp_A0 = (A_v_miss{1,vv}+A_v_miss{1,vv}')/2;
    if isnormalization == 0
        temp_D0 = diag(sum(temp_A0));
        temp_L0 = temp_D0 - temp_A0;
        [temp_F_v,~,~] = svds(temp_L0,num_class);
    else
        temp_D10 = diag(1./sqrt(sum(temp_A0, 2)+eps));
        A_v_miss{1,vv} = temp_D10*temp_A0*temp_D10;
        [temp_F_v,~,~] = svds(A_v_miss{1,vv},num_class);
    end
    
    F_v{1,vv} = temp_F_v;
    A_v{1,vv} = A_v_miss{1,vv};
end

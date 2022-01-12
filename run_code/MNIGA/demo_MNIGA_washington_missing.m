% clc
clear

addpath(genpath('../../generate_dataset'))
addpath(genpath('../../imnga_codes'))


%% load data 
load('washington_missing.mat')
clear X X_id Y
X = X_missing07;
X_id = X_missing07_id;
Y = Y_missing07;

%% Preprocessing and parameter setting
max_iter_out = 50;
paras_space = [10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4];
lambda1 = paras_space(4); lambda2 = paras_space(5); lambda3 = paras_space(9);
alpha = 1;
%% 
num_view = size(X,2);
for group_id = 1:15  % 15 randomly generated incomplete groups
    for vv = 1:num_view
        % X{group_id,vv} = whiteningPCA(X{group_id,vv}, 25);
        if size(X{group_id,vv},2)>250
            X{group_id,vv} = NormalizeFea(X{group_id,vv},1);
        end
        X_miss{1,vv} = X{group_id,vv};
        X_miss_id{1,vv} = X_id{group_id,vv};
    end
    Y_miss = Y{group_id};

    %% Initialization
    num_view = size(X_miss,2);
    num_class = max(Y_miss);
    num_samples = size(Y_miss,1);
    
    [A_v, F_v, G_v, M_v, num_samples_v] = ...
        initialization_without_compensate(X_miss,X_miss_id,num_class,num_samples);
   
    %% run IMNGA 
    [H, ~, ~, A, ~, result_imcan1_all, delte_v_all,Tobj] =  ...
            Imc_IMNGA(X_miss, Y_miss, A_v, F_v, G_v, lambda1, ...
            lambda2, lambda3, max_iter_out, num_samples_v, alpha);

    [non ,Y_Hpred] = max(H,[],2);
    % [ACC MI hat Purity]
    result_imcan1_paras(group_id,:) = ClusteringMeasure(Y_miss, Y_Hpred);
    result_imcan1_paras2(group_id,:) = max(result_imcan1_all,[],1);
end
[mean(result_imcan1_paras2); std(result_imcan1_paras2)]'
[mean(result_imcan1_paras); std(result_imcan1_paras)]';


%% selest parameters
% num_view = size(X,2);
% for group_id = 1:15  % 15 randomly generated incomplete groups
%     for vv = 1:num_view
%         % X{group_id,vv} = whiteningPCA(X{group_id,vv}, 25);
%         if size(X{group_id,vv},2)>250
%             X{group_id,vv} = NormalizeFea(X{group_id,vv},1);
%         end
%         X_miss{1,vv} = X{group_id,vv};
%         X_miss_id{1,vv} = X_id{group_id,vv};
%     end
%     Y_miss = Y{group_id};
% 
%     %% Initialization
%     num_view = size(X_miss,2);
%     num_class = max(Y_miss);
%     num_samples = size(Y_miss,1);
%     
%     [A_v, F_v, G_v, M_v, num_samples_v] = ...
%         initialization_without_compensate(X_miss,X_miss_id,num_class,num_samples);
%    
%     %% run IMNGA 
%     paras_space = [10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2 10^3 10^4];
%     for ii = 1:5  % lambda_2
%         for jj = 1:5   % lambda_1
%             for kk = 4:9
%                 [group_id ii jj kk]
%                 [H, ~, ~, A, ~, result_imcan1_all, delte_v_all,Tobj] =  ...
%                     Imc_IMNGA(X_miss, Y_miss, A_v, F_v, G_v, paras_space(ii), ...
%                     paras_space(jj), paras_space(kk), max_iter_out, num_samples_v, alpha);
%                 [non ,Y_Hpred] = max(H,[],2);
%                 % [ACC MI hat Purity]
%                 result_imcan1_paras(group_id,ii,jj,kk,:) = ClusteringMeasure(Y_miss, Y_Hpred);
%                 result_imcan1_paras2(group_id,ii,jj,kk,:) = max(result_imcan1_all,[],1);
%             end
%         end
%     end
% end
% 
% size(result_imcan1_paras2,1)
% results1_add = zeros(size(result_imcan1_paras2,1),5,size(result_imcan1_paras2,3),9);
% for kk =1:3
%     results1_add = results1_add + result_imcan1_paras2(:,:,:,:,kk);
% end
% results1_best1 = zeros(size(results1_add,1),3);
% results1_sum = [];
% results1_sum(:,:,:) = sum(results1_add,1);
% for ii=1:5
%     for jj=1:size(result_imcan1_paras2,3)
%         for kk=1:9
%             if results1_sum(ii,jj,kk)==max(max(max(max(max(results1_sum)))))
%                 [ii,jj,kk]
%                 index_x1=ii;index_x2=jj;index_x3=kk;
%             end
%         end
%     end
% end
% results1_best1(:,:) = result_imcan1_paras2(:,index_x1,index_x2,index_x3,:);
% [mean(results1_best1,1); std(results1_best1,1)]'
% results1_add = zeros(size(results1_add,1),5,size(result_imcan1_paras2,3),9);
% for kk =1:3
%     results1_add = results1_add + result_imcan1_paras(:,:,:,:,kk);
% end
% results1_best2 = zeros(size(results1_add,1),3);
% [~,index] = max(results1_add,[],2);
% for iii=1:size(results1_add,1)
%     results1_sum(:,:,:,:) = results1_add(iii,:,:,:,:);
%     for ii=1:5
%         for jj=1:size(result_imcan1_paras2,3)
%             for kk=1:9
%                 if results1_sum(ii,jj,kk)==max(max(max(max(max(results1_sum)))))
%                     index_x1=ii;index_x2=jj;index_x3=kk;
%                     [iii ii,jj,kk];
%                     results1_sum(ii,jj,kk);
%                 end
%             end
%         end
%     end
%     results1_best2(iii,:) = result_imcan1_paras(iii,index_x1,index_x2,index_x3,:);
% end
% [mean(results1_best2,1); std(results1_best2,1)]';
% results1_best2;
% [mean(results1_best2,1); std(results1_best2,1);mean(results1_best1,1); std(results1_best1,1);]'

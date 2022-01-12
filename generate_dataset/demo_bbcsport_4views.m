clear
class = load('bbcsport.clist');
class = sort(class);
Y = floor(class'/10000);
for i=1:max(Y)
    num_class(i) = sum(Y(:)==i);
end
num_class0 = [0 num_class];
view1 = load('bbcsport_seg1of4.mtx');
view2 = load('bbcsport_seg2of4.mtx');
view3 = load('bbcsport_seg3of4.mtx');
view4 = load('bbcsport_seg4of4.mtx');
for ii=1:size(view1)
    X1(view1(ii,1),view1(ii,2))=view1(ii,3);
end
X1 = X1';
for ii=1:size(view2)
    X2(view2(ii,1),view2(ii,2))=view2(ii,3);
end
X2 = X2';
for ii=1:size(view3)
    X3(view3(ii,1),view3(ii,2))=view3(ii,3);
end
X3 = X3';
for ii=1:size(view4)
    X4(view4(ii,1),view4(ii,2))=view4(ii,3);
end
X4 = X4';

X1_missing = load('bbcsport_seg1of4.txt');
X1_missing_class = floor(X1_missing/10000);
X1_missing_rem = rem(X1_missing,10000);
for ii=1:length(X1_missing)
    X1_missing(ii) = X1_missing_rem(ii) + sum(num_class0(1:X1_missing_class(ii)));
end
X2_missing = load('bbcsport_seg2of4.txt');
X2_missing_class = floor(X2_missing/10000);
X2_missing_rem = rem(X2_missing,10000);
for ii=1:length(X2_missing)
    X2_missing(ii) = X2_missing_rem(ii) + sum(num_class0(1:X2_missing_class(ii)));
end
X3_missing = load('bbcsport_seg3of4.txt');
X3_missing_class = floor(X3_missing/10000);
X3_missing_rem = rem(X3_missing,10000);
for ii=1:length(X3_missing)
    X3_missing(ii) = X3_missing_rem(ii) + sum(num_class0(1:X3_missing_class(ii)));
end
X4_missing = load('bbcsport_seg4of4.txt');
X4_missing_class = floor(X4_missing/10000);
X4_missing_rem = rem(X4_missing,10000);
for ii=1:length(X4_missing)
    X4_missing(ii) = X4_missing_rem(ii) + sum(num_class0(1:X4_missing_class(ii)));
end
X{1,1} = X1; X{1,2} = X2; X{1,3} = X3; X{1,4} = X4;
X_id{1,1} = X1_missing; X_id{1,2} = X2_missing;
X_id{1,3} = X3_missing; X_id{1,4} = X4_missing;

temp1=[X_id{1,1}; X_id{1,2}; X_id{1,3}; X_id{1,4}];
temp2 = tabulate(temp1(:));
temp3 = (temp2(:,2)==4);
num_samples = size(Y,1);
X_new=[];
for vv=1:4
    temp_X_v_id = X_id{1,vv};
    num_samples_v(vv,1) = size(temp_X_v_id,1);
    G_v{1,vv} = zeros(num_samples_v(vv,1),num_samples);
    for sample_ii = 1:num_samples_v(vv,1)
        G_v{1,vv}(sample_ii,temp_X_v_id(sample_ii))=1;
    end
    G2_v{1,vv} = G_v{1,vv}(:,temp3);
    X_new{1,vv} = G2_v{1,vv}' * X{1,vv};
    X_id_new{1,vv} = [1:116]';
end
Y_new = Y(temp3,1);
X = X_new;
X_id = X_id_new;
Y = Y_new;
clear G2_v G_v temp1 temp2 temp3 temp_X_v_id vv num_samples X_new X_id_new Y_new num_samples_v sample_ii


ii = 1;
while ii<16
    contact01 = []; contact03 = []; contact05 = [];
    for jj=1:4
        length_ = size(X_id{1,jj},1);
        temp_rand = randperm(length_);
        X_missing01{ii,jj} = X{1,jj}(sort(temp_rand(1:ceil(length_*0.9))),:);
        X_missing01_id{ii,jj} = X_id{1,jj}(sort(temp_rand(1:ceil(length_*0.9))),:);
        contact01 = [contact01; X_missing01_id{ii,jj}(:)];
        X_missing03{ii,jj} = X{1,jj}(sort(temp_rand(1:ceil(length_*0.7))),:);
        X_missing03_id{ii,jj} = X_id{1,jj}(sort(temp_rand(1:ceil(length_*0.7))),:);
        contact03 = [contact03; X_missing03_id{ii,jj}(:)];
        X_missing05{ii,jj} = X{1,jj}(sort(temp_rand(1:ceil(length_*0.5))),:);
        X_missing05_id{ii,jj} = X_id{1,jj}(sort(temp_rand(1:ceil(length_*0.5))),:);
        contact05 = [contact05; X_missing05_id{ii,jj}(:)];
    end
    
    if length(unique(contact01))==length(Y)
        Y_missing01{ii,1} = Y;
    else
        temp_id = sort(unique(contact01));
        Y_missing01{ii,1} = Y(temp_id,1);
        for kk = 1:length(temp_id)
            temp_id_inv(temp_id(kk)) =  kk;
        end
        for jj=1:4
            for kk = 1:length(X_missing01_id{ii,jj})
                X_missing01_id{ii,jj}(kk,1) = temp_id_inv(X_missing01_id{ii,jj}(kk,1));
            end
        end
    end
    
    if length(unique(contact03))==length(Y)
        Y_missing03{ii,1} = Y;
    else
        temp_id = sort(unique(contact03));
        Y_missing03{ii,1} = Y(temp_id,1);
        for kk = 1:length(temp_id)
            temp_id_inv(temp_id(kk)) =  kk;
        end
        for jj=1:4
            for kk = 1:length(X_missing03_id{ii,jj})
                X_missing03_id{ii,jj}(kk,1) = temp_id_inv(X_missing03_id{ii,jj}(kk,1));
            end
        end
    end
    
    if length(unique(contact05))==length(Y)
        Y_missing05{ii,1} = Y;
    else
        temp_id = sort(unique(contact05));
        Y_missing05{ii,1} = Y(temp_id,1);
        for kk = 1:length(temp_id)
            temp_id_inv(temp_id(kk)) =  kk;
        end
        for jj=1:4
            for kk = 1:length(X_missing05_id{ii,jj})
                X_missing05_id{ii,jj}(kk,1) = temp_id_inv(X_missing05_id{ii,jj}(kk,1));
            end
        end
    end
    ii = ii + 1;
end

save BBCsport_4views116_missing X_id X Y X_missing01 X_missing01_id Y_missing01 ...
    X_missing03 X_missing03_id Y_missing03 X_missing05 X_missing05_id Y_missing05

%%
tempY = Y_missing03;
for ii=1:15
    ii
    length(tempY{ii,1})
    aaa = tabulate(tempY{ii,1}(:));
    aaa(:,2)
end

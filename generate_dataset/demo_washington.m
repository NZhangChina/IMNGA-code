% clear
view1 = load('washington_cites.mtx');
view2 = load('washington_content.mtx');
for ii=2:size(view1)
    X1(view1(ii,1),view1(ii,2))=view1(ii,3);
end
X{1,1} = X1;
for ii=2:size(view2)
    X2(view2(ii,1),view2(ii,2))=view2(ii,3);
end
X{1,2} = X2;

num_view = 2;
Y = load('washington_act.txt');

temp_id = 1:length(Y);
X_id{1,1} = temp_id'; X_id{1,2} = temp_id';
%%
ii = 1;
while ii<16
    contact03 = []; contact05 = []; contact07 = []; contact09 = [];
    length_ = size(X_id{1,1},1);
    temp_rand = randperm(length_);
    for jj=1:num_view
        rand_id03{1,jj} = [temp_rand(1:ceil(length_*0.3)) temp_rand(ceil(length_*0.3+length_*(jj-1)*0.7/num_view)+1:ceil(length_*0.3+length_*jj*0.7/num_view))];
        rand_id05{1,jj} = [temp_rand(1:ceil(length_*0.5)) temp_rand(ceil(length_*0.5+length_*(jj-1)*0.5/num_view)+1:ceil(length_*0.5+length_*jj*0.5/num_view))];
        rand_id07{1,jj} = [temp_rand(1:ceil(length_*0.7)) temp_rand(ceil(length_*0.7+length_*(jj-1)*0.3/num_view)+1:ceil(length_*0.7+length_*jj*0.3/num_view))];
        rand_id09{1,jj} = [temp_rand(1:ceil(length_*0.9)) temp_rand(ceil(length_*0.9+length_*(jj-1)*0.1/num_view)+1:ceil(length_*0.9+length_*jj*0.1/num_view))];
    end
    for jj=1:num_view
        X_missing03{ii,jj} = X{1,jj}(sort(rand_id03{1,jj}),:);
        X_missing03_id{ii,jj} = X_id{1,jj}(sort(rand_id03{1,jj}),:);
        
        X_missing05{ii,jj} = X{1,jj}(sort(rand_id05{1,jj}),:);
        X_missing05_id{ii,jj} = X_id{1,jj}(sort(rand_id05{1,jj}),:);
        
        X_missing07{ii,jj} = X{1,jj}(sort(rand_id07{1,jj}),:);
        X_missing07_id{ii,jj} = X_id{1,jj}(sort(rand_id07{1,jj}),:);
        
        X_missing09{ii,jj} = X{1,jj}(sort(rand_id09{1,jj}),:);
        X_missing09_id{ii,jj} = X_id{1,jj}(sort(rand_id09{1,jj}),:);

    end
    Y_missing03{ii,1} = Y;
    Y_missing05{ii,1} = Y;
    Y_missing07{ii,1} = Y;
    Y_missing09{ii,1} = Y;
    ii = ii + 1;
end

save washington_missing X_id X Y X_missing03 X_missing03_id Y_missing03 ...
      X_missing05 X_missing05_id Y_missing05 X_missing07 X_missing07_id ...
      Y_missing07 X_missing09 X_missing09_id Y_missing09



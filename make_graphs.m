%hist(scores)
not_shared_target_indices = find(indicators == 0);
shared_target_indices = find(indicators == 1);
scores_1 = scores(shared_target_indices,:);
scores_0 = scores(not_shared_target_indices,:);
%hist(scores_0);
hist(scores_1);

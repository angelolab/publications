function [mean_peak_loc] = peak_loc(trajectories)
[~,I] = max(trajectories,[],2);
mean_peak_loc = mean(I);
end


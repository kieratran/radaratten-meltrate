function [Na_pts, m_pts] = off_bound_detection(attenuation, melt)
% Finding the cut-off bounds where attenuation is no longer reflecting
% changes in ice temperature, hence basal melt rates

second_der = diff(attenuation, 2); % second derivative of attenuation rates
x = melt(3:end); % matching basal melt rate matrix to second derivative matrix
y = attenuation(3:end); % matching attenuation rate matrix to second derivative matrix
tol = 0.05; % flateau threshold
nonflat_idx = find(abs(second_der) > tol); % all index values above the threshold
if ~isempty(nonflat_idx)
    Na_pts(2,:) = y(nonflat_idx(end)); % lower bound of attenuation rate
    m_pts(2,:) = x(nonflat_idx(end)); % upper bound of basal melt rate
    Na_pts(1,:) = y(nonflat_idx(1)); % upper bound of attenuation rate
    m_pts(1,:) = x(nonflat_idx(1)); % lower bound of basal melt rate
end

end
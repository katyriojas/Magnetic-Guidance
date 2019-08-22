function [T, FREmatch] = point_register_without_correspondence(X,Y,plot_flag)
% Outputs transformation from matrix of X pnts to matrix of Y points
% Does not require point correspondance

% Each column of X and each column of Y represents a
% K-dimensional point, where K is the number of rows.
% X must have more or the same number of points (columns) as Y.

%%
if nargin > 3
    error('Error. Too many inputs.')
elseif nargin < 2
    error('Error. Need two sets of points as inputs.')
elseif nargin < 3
    plot_flag = 0; 
end

if (plot_flag==1) && (size(Y,1)~=3)
    error('Error. Can only plot 3D points.');
end

%% Brute force registration
% create matrix where each row is a unique ordered set
combinations = nchoosek(1:size(X,2), size(Y,2));
permutations = perms(combinations(1,:)); % initialize permutations matrix
for i=2:size(combinations,1)
    permutations = [permutations; perms(combinations(i,:))];
end

% try all possible registrations
FRE = zeros(size(permutations,1),1);
for i=1:size(permutations,1)
    [~,~,FRE(i)] = point_register(X(:,permutations(i,:)), Y);
end

%% find the transformation for the corresponding set
[FREmin,index] = min(FRE); % lowest FRE is the correct set
XMatched = X(:, permutations(index,:));
[Rmatch, tmatch, FREmatch] = point_register(XMatched, Y);
T = [Rmatch, tmatch; [0 0 0 1]];

%% plot
if plot_flag
    figure
    hold on
    
    % Put X points in Y frame
    XRegistered = Rmatch*X + tmatch*ones(1,size(X,2));
    XMatchedRegistered = Rmatch*XMatched + tmatch*ones(1,size(Y,2));
    
    % plot
    drawPolyline3d(Y','closed','-g','LineWidth',1.3)
    drawPolyline3d(XMatchedRegistered','closed','--r','LineWidth',1.3)
    drawPoint3d(XRegistered','ob')
    
    legend('Y', 'XMatched', 'X')
    result = sprintf('Best FRE = %1.3g', FREmin);
    title(result)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    view(45,45)
    axis equal vis3d
    
    figure
    FREsort = sort(FRE);
    plot(FREsort);
    hold on
    xlabel('permutation')
    ylabel('FRE')
    plot(1,FREsort(1)','og');
    plot(2,FREsort(2)','*m');
    title(sprintf('Best FRE = %1.3g, Next = %1.3g',FREsort(1),FREsort(2)))
end

end

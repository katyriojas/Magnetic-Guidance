function [R,t,FRE,FREcomponents] = point_register(X,Y,w,n_t)
% [R,t,FRE,FREcomponents] = POINT_REGISTER(X,Y,w,n_t)
% Find R,t such that R*X + t*ones(1,N) ~ Y,
% where N = number of columns of X and Y. Optional weightings may be
% specified in the vector w. n_t is the number of points which contribute
% to the calculation of the translation vector t. If n_t is zero,
% then t is all zeros. If n_t is omitted, then all points contribute to t.
%
% Author: J. Michael Fitzpatrick
% Creation: Fall 2001 (for my CS 357 course, Image Processing)
%
% Each column of X and each column of Y represents a
% K-dimensional point, where K is the number of rows.
% X and Y must have the same number of points (columns).
% FRE = RMS fiducial registration error.
% FREcomponents = R*X + t*ones(1,N) - Y.
% Uses algorithm 8.1 and notation of pp. 469-70 of
% J. Michael Fitzpatrick, Derek L. G. Hill, and Calvin R. Maurer, Jr.,
% "Image Registration", Chapter 8 of "Handbook of Medical Imaging,
% Volume 2, Medical Image Processing and Analysis",
% Milan Sonka and J. Michael Fitzpatrick, eds., SPIE Press,
% Bellingham, Wa, 2000.
%
% Sep 18, 2005:
% Modified by JMF to allow arb. value for K. In
% particular, it is able to handle K = 1 or 2, as well as K = 3, but it is
% written to handle any number of dimensions. The method of
% insuring that R is proper (by means of the diagonal matrix D) has
% been proven to find the optimum transformation only up to 3 dimensions.
%
% May 8, 2007:
% Modified by JMF to allow N < K. This is a trivial
% modification: Removal of a check for N<K. The result is,
% for example, that a registration can be found for two points
% in three-dimensional space (N = 2, K = 3).
%
% May 8, 2007:
% Also modified by JMF to add the optional fourth argument.
% This argument is the number of points (beginning with the first)
% for which translation is to be calculated and used. For the rest
% of the points only rotation is used. This feature makes new
% applications possible. For example, if it is desired to use
% a set of points for which the first subset are real points,
% but the second subset are unit vectors representing direction
% only, then n_t is the number of real points.
%
% February 18, 2010
% Modified by JMF to insert the optional third argument, specifying
% (isotropic) point weightings. The algorithm is the same one referenced
% above. This option was tested by Ramya Balachandran on an earlier version
% of point_register, which allowed a set of squared weights to be given,
% but in this version, the third argument is a vector of unsquared weights.
% These weights are normalized so that their squares sum to one. This means
% that, for example, if each weight is multiplied by 100, FRE will not
% change.
%
% September 9, 2011:
% Modified by JMF to correct a bug in the calculation of the translation
% vector for the case of n_t < N. 

% Validate and preprocess input:
if nargin < 2
   error('At least two input arguments are required.');
end
[K N] = size(X);
[K_Y N_Y] = size(Y);
if K ~= K_Y
   error('X and Y must have the same number of rows.')
elseif (K ~= K_Y) || (N ~= N_Y)
   error('X and Y must have the same number of columns.');
end
if nargin >= 3
   if ~isvector(w)
      error('w must be a vector');
   end
   if length(w) ~= N
      error('length of w must equal number of columns in X and Y');
   end
   w = w(:); % force it to be a column vector
else
   w = ones(N,1);
end
if nargin < 4
   n_t = N; % User wants to use all points to compute translation
elseif n_t > N
   n_t = N; % because larger numbers make no sense
elseif n_t < 0
   n_t = 0;  % because negative numbers make no sense
end

% Calculate normalized weights and centroids:
w_sqr = w.^2;
sum_w_sqr = sum(w_sqr);
w_sqr_normed = w_sqr/sum_w_sqr;
if n_t == 0
   Xbar = zeros(K,1); Ybar = zeros(K,1); % dummy centroids
else
   w_sqr_n_t = w_sqr(1:n_t);
   w_sqr_n_t_normed = w_sqr_n_t/sum(w_sqr_n_t); % normed over first n_t weights
   Xbar = X(:,1:n_t)*w_sqr_n_t_normed;  % X weighted centroid of first n_t points
   Ybar = Y(:,1:n_t)*w_sqr_n_t_normed;  % Y weighted centroid of first n_t points
end

% Demean the first n_t points:
Xtilde = X - [repmat(Xbar,1,n_t),zeros(K,N-n_t)];
Ytilde = Y - [repmat(Ybar,1,n_t),zeros(K,N-n_t)];

% Calculate the rotation:
H = Xtilde*diag(w_sqr_normed)*Ytilde';  % cross covariance matrix
[U S V] = svd(H);    % U*S*V' = H (change S to ~ for later Matlab versions)
D = diag([ones(1,K-1), det(V*U)]); % used to insure that R is proper
R = V*D*U';

% Calculate the translation:
t = Ybar-R*Xbar;

% Finish generating output:
if nargout >= 3
   if n_t ~= N
      T = [repmat(t,1,n_t),zeros(K,N-n_t)];
   else
      T = repmat(t,1,N);
   end
   FREcomponents = repmat(w',K,1).*(R*X + T - Y);
   FRE = sqrt(FREcomponents(:)'*FREcomponents(:)/sum_w_sqr);
end






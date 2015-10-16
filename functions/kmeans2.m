function [idx, C] = kmeans2(X)
%KMEANS K-means clustering.
%   IDX = KMEANS(X, 2) partitions the points in the N-by-P data matrix X
%   into K clusters.  This partition minimizes the sum, over all clusters, of
%   the within-cluster sums of point-to-cluster-centroid distances.  Rows of X
%   correspond to points, columns correspond to variables.  Note: when X is a
%   vector, KMEANS treats it as an N-by-1 data matrix, regardless of its
%   orientation.  KMEANS returns an N-by-1 vector IDX containing the cluster
%   indices of each point.  By default, KMEANS uses squared Euclidean
%   distances.
%
%   KMEANS treats NaNs as missing data, and ignores any rows of X that
%   contain NaNs.
%
%   [IDX, C] = KMEANS(X, K) returns the K cluster centroid locations in
%   the K-by-P matrix C.
%
%   Example:
%
%       X = [randn(20,2)+ones(20,2); randn(20,2)-ones(20,2)];
%       [cidx, ctrs] = kmeans(X, 2);
%       plot(X(cidx==1,1),X(cidx==1,2),'r.', ...
%            X(cidx==2,1),X(cidx==2,2),'b.', ctrs(:,1),ctrs(:,2),'kx');
%
%   See also LINKAGE, CLUSTERDATA, SILHOUETTE.

%   KMEANS uses a two-phase iterative algorithm to minimize the sum of
%   point-to-centroid distances, summed over all K clusters.  The first phase
%   uses what the literature often describes as "batch" updates, where each
%   iteration consists of reassigning points to their nearest cluster
%   centroid, all at once, followed by recalculation of cluster centroids.
%   This phase occasionally (especially for small data sets) does not converge
%   to solution that is a local minimum, i.e., a partition of the data where
%   moving any single point to a different cluster increases the total sum of
%   distances.  Thus, the batch phase be thought of as providing a fast but
%   potentially only approximate solution as a starting point for the second
%   phase.  The second phase uses what the literature often describes as
%   "on-line" updates, where points are individually reassigned if doing so
%   will reduce the sum of distances, and cluster centroids are recomputed
%   after each reassignment.  Each iteration during this second phase consists
%   of one pass though all the points.  The on-line phase will converge to a
%   local minimum, although there may be other local minima with lower total
%   sum of distances.  The problem of finding the global minimum can only be
%   solved in general by an exhaustive (or clever, or lucky) choice of
%   starting points, but using several replicates with random starting points
%   typically results in a solution that is a global minimum.
%
% References:
%
%   [1] Seber, G.A.F. (1984) Multivariate Observations, Wiley, New York.
%   [2] Spath, H. (1985) Cluster Dissection and Analysis: Theory, FORTRAN
%       Programs, Examples, translated by J. Goldschmidt, Halsted Press,
%       New York.

%   Copyright 1993-2014 The MathWorks, Inc.

% n points in p dimensional space
n = size(X,1);

maxit = 100;


% case {'plus','kmeans++'}
% Select the first seed by sampling uniformly at random
C(1,:) = X(ceil(end*rand),:);

% Select the rest of the seeds by a probabilistic model
sampleProbability = (X(:,1) - C(1)).^2 + (X(:,2) - C(2)).^2;
denominator = sum(sampleProbability);
sampleProbability = sampleProbability/denominator;
edges = min([0; cumsum(sampleProbability)],1); % protect against accumulated round-off
edges(end) = 1; % get the upper edge exact
ps = rand;
C(2,:) = X(edges(1:end-1)<=ps&ps<edges(2:end),:);

% Compute the distance from every point to each cluster centroid and the
% initial assignment of points to clusters
D = distfun(X, C);
[~, idx] = min(D, [], 2);
[C, m] = gcentroids(X, idx);

% Begin phase one:  batch reassignments
%------------------------------------------------------------------
% Every point moved, every cluster will need an update
prevtotsumD = Inf;
iter = 0;
while true
    iter = iter + 1;
    
    % Calculate the new cluster centroids and counts, and update the
    % distance from every point to those new cluster centroids
    Clast = C;
    mlast = m;
    D = distfun(X, C);
    
    % Deal with clusters that have just lost all their members
    if any(m==0)
        i = find(m==0);
        d = D((idx-1)*n + (1:n)'); % use newly updated distances
        
        % Find the point furthest away from its current cluster.
        % Take that point out of its cluster and use it to create
        % a new singleton cluster to replace the empty one.
        [~, lonely] = max(d);
        from = idx(lonely); % taking from this cluster
        if m(from) < 2
            % In the very unusual event that the cluster had only
            % one member, pick any other non-singleton point.
            from = find(m>1,1,'first');
            lonely = find(idx==from,1,'first');
        end
        idx(lonely) = i;
        
        % Update clusters from which points are taken
        [C, m] = gcentroids(X, idx);
        D = distfun(X, C);
    end
    
    % Compute the total sum of distances for the current configuration.
    totsumD = sum(D((idx-1)*n + (1:n)'));
    % Test for a cycle: if objective is not decreased, back out
    % the last step and move on to the single update phase
    if prevtotsumD <= totsumD
        idx = previdx;
        C = Clast;
        m = mlast;
        iter = iter - 1;
        break;
    end
    if iter >= maxit
        break;
    end
    
    % Determine closest cluster for each point and reassign points to clusters
    previdx = idx;
    prevtotsumD = totsumD;
    nidx = (D(:,2)<D(:,1))+1;
    
    % Determine which points moved
    moved = nidx ~= previdx;
    if any(moved)
        % Resolve ties in favor of not moving
        moved(D(moved,1)==D(moved,2)) = false;
    end
    if ~any(moved)
        break;
    end
    idx(moved) = nidx(moved);
    [C, m] = updateCenters(C,X,m,nidx,moved);
    
    % DN: don't update changed. Since we are only looking for two clusters,
    % if any point changed, both clusters have to be updated.
    % Find clusters that gained or lost members
    % changed = unique([idx(moved); previdx(moved)])';
    
end % phase one
%------------------------------------------------------------------

% Begin phase two:  single reassignments
%------------------------------------------------------------------
lastmoved = 0;
nummoved = 0;
converged = false;
while iter < maxit
    % Calculate distances to each cluster from each point, and the
    % potential change in total sum of errors for adding or removing
    % each point from each cluster.  Clusters that have not changed
    % membership need not be updated.
    %
    % Singleton clusters are a special case for the sum of dists
    % calculation.  Removing their only point is never best, so the
    % reassignment criterion had better guarantee that a singleton
    % point will stay in its own cluster.  Happily, we get
    % Del(i,idx(i)) == 0 automatically for them.
    Del = distfun(X, C);
    mbrs = idx==1;
    sgn = 1 - 2*mbrs; % -1 for members, 1 for nonmembers
    if m(1) == 1
        sgn(mbrs) = 0; % prevent divide-by-zero for singleton mbrs
    end
    Del(:,1) = (m(1) ./ (m(1) + sgn)) .* Del(:,1);
    % same for cluster 2
    sgn = -sgn; % -1 for members, 1 for nonmembers
    if m(2) == 1
        sgn(~mbrs) = 0; % prevent divide-by-zero for singleton mbrs
    end
    Del(:,2) = (m(2) ./ (m(2) + sgn)) .* Del(:,2);
    
    % Determine best possible move, if any, for each point.  Next we
    % will pick one from those that actually did move.
    previdx = idx;
    nidx = (Del(:,2)<Del(:,1))+1;
    moved = find(previdx ~= nidx);
    if ~isempty(moved)
        % Resolve ties in favor of not moving
        moved(Del(moved,1)==Del(moved,2)) = [];
    end
    if isempty(moved)
        converged = true;
        break;
    end
    
    % Pick the next move in cyclic order
    moved = mod(min(mod(moved - lastmoved - 1, n) + lastmoved), n) + 1;
    
    % If we've gone once through all the points, that's an iteration
    if moved <= lastmoved
        iter = iter + 1;
        if iter >= maxit, break; end
        nummoved = 0;
    end
    nummoved = nummoved + 1;
    lastmoved = moved;
    
    oidx = idx(moved);
    nidx = nidx(moved);
    totsumD = totsumD + Del(moved,nidx) - Del(moved,oidx);
    
    % Update the cluster index vector, and the old and new cluster
    % counts and centroids
    idx(moved) = nidx;
    m(nidx) = m(nidx) + 1;
    m(oidx) = m(oidx) - 1;
    C(nidx,:) = C(nidx,:) + (X(moved,:) - C(nidx,:)) / m(nidx);
    C(oidx,:) = C(oidx,:) - (X(moved,:) - C(oidx,:)) / m(oidx);
end % phase two
%------------------------------------------------------------------
if ~converged
    warning(message('stats:kmeans:FailedToConverge', maxit));
end





%------------------------------------------------------------------
function D = distfun(X, C)
%DISTFUN Calculate point to cluster centroid distances.
D = [(X(:,1) - C(1,1)).^2 + (X(:,2) - C(1,2)).^2,  (X(:,1) - C(2,1)).^2 + (X(:,2) - C(2,2)).^2];

%------------------------------------------------------------------
function [centroids, counts] = gcentroids(X, index)
%GCENTROIDS Centroids and counts stratified by group.
n = size(X,1);
centroids = NaN(2);
counts = zeros(2,1);

members = index == 1;
counts(1) = sum(members);
centroids(1,:) = sum(X(members,:),1) / counts(1);
counts(2) = n-counts(1);
centroids(2,:) = sum(X(~members,:),1) / counts(2);

%------------------------------------------------------------------
function [C, m] = updateCenters(C,X,m,nidx,moved)
idx    = find(moved);
nidx   = nidx(moved);
oidx   = 3-nidx;  %nidx==1 -> oidx=2, nidx==2 -> oidx=1
nmoved = length(idx);
for p=1:nmoved
    m(nidx(p)) = m(nidx(p))+1;
    m(oidx(p)) = m(oidx(p))-1;
    C(nidx(p),:) = C(nidx(p),:) + (X(idx(p),:)-C(nidx(p),:))/m(nidx(p));
    C(oidx(p),:) = C(oidx(p),:) - (X(idx(p),:)-C(oidx(p),:))/m(oidx(p));
end
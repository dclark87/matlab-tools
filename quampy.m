%% QUAMPY Quanternion multiplication function
% PQ = QUAMPY(P,Q) returns the quaternion multiplication of two quaternions
% or if either is just a 3D vector, it will append a zero the beginning of the
% vector to form a pure quaternion in the form of p = [0 p(1) p(2) p(3)]
%
% where p = [p0 pi pj pk]; q = [q0 qi qj qk]
% and
% pq = p*q where * is the quaternion multiplication
%    = p0q0 - <p,q> + p0q + q0p + p x q
%
% Author: Daniel Clark, 2013
%%
function pq = quampy(p,q)
% Check vector size, p
if numel(p) ~= 4 && numel(p) ~= 3
    error('Number of elements in 1st vector needs to be 3 or 4');
else
    p = reshape(p,1,numel(p));
end
if numel(p) == 3, p = [0 p]; end
% Check vector size, q
if numel(q) ~= 4 && numel(q) ~= 3
    error('Number of elements in 2nd vector needs to be 3 or 4');
else
    q = reshape(q,1,numel(q));
end
if numel(q) == 3, q = [0 q]; end

pq0 = p(1)*q(1) - dot(p(2:end),q(2:end));
pq = p(1)*q(2:end) + q(1)*p(2:end) + cross(p(2:end),q(2:end));
pq = [pq0 pq];
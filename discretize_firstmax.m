%
% Author: Jan Lellmann, j.lellmann@damtp.cam.ac.uk
%

function result = discretize_firstmax(uinput, dims, gamma)
%DISCRETIZE_FIRSTMAX Discretizes (->unit vectors) u by
%   assigning to each vector the first unit vector belonging
%   to a maximal component. Works over the components of u.
%   If the (dims.components)-dimensional vector gamma is given, each
%   vector in u is scaled elementwise by gamma before taking the max.
%   There are no assumptions on gamma, but is reasonable to force the
%   elements to be non-negative and sum to one (as adding the same constant
%   to all elements in gamma will not change the result).

    u = reshape(uinput, [dims.nimage dims.components]);
    
    if (nargin > 2)
        u = u .* repmat(gamma(:)',dims.nimage,1);
    end

    [tmp,idx] = max(u,[],2); % indices of maximal elements along 2nd dimension
    result = zeros(dims.nimage, dims.components);
    result((1:dims.nimage)' + (idx - 1)*dims.nimage) = 1; % set i-th row to idx(i)-th unit vector    
    result = reshape(result, size(uinput));    
end

% previous code (slow, not using index variant of max)
%    result = zeros(size(u));   
%    for i = 1:size(u,1) %SLOW
%        idx = find(u(i,:) == max(u(i,:)));
%        result(i,idx(1)) = 1;
%    end

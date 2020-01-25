%
% Author: Jan Lellmann, j.lellmann@damtp.cam.ac.uk
%

% Represents the data term <u,s> for fixed s.
% The Lipschitz constant is zero (valid for any norm!).
function result = dataterm_linear(s,dims)
    if (nargin < 2)
        dims = dimensions(numel(s),1);
    end
    result = [];
    result.evaluate = @(u)(dataterm_linear_evaluate(u,s,dims));
    result.gradient = @(u)(s(:)');
    result.gradient_lipschitz = 0;
    result.islinear = true;    
	
	result.dims = dims; % this is private and should not be used
end

function [result,resultlocal] = dataterm_linear_evaluate(u,s,dims)
    r = sum(reshape(u(:) .* s(:),dims.nimage,dims.components),2);
    result = sum(r);
    if (nargout > 1)
        resultlocal = r;
    end
end

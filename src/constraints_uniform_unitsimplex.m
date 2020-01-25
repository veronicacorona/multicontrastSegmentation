%
% Author: Jan Lellmann, j.lellmann@damtp.cam.ac.uk
%

% If one argument is given, projects each vector separately onto the unit
% simplex.

% (obsolete) If basis is supplied, it must be a (k x dims.components)
% matrix; in this case each vector is projected on the intersection of the
% unit simplex with the orthogonal complement of the kernel of basis.
% Set basis = [] to disable.

% mask and border can be used to specify "ghost" cells that are fixed.
% mask ist an index vector of size dims.ntotal x 1; border is a
% numerical t x 1 vector with t the number of "true" elements in mask.
% At each iteration, these cells are reset via u(mask) = border(:).
% If mask contains an entry for one component at a pixel, it must contain
% entries for all other components at that pixel as well; otherwise the
% behavior is undefined!
function result = constraints_uniform_unitsimplex(dims,basis,mask,border)
    if (nargin > 1 && ~isempty(basis))
        error('Not implemented yet!');
    end
    if (nargin >= 2 && nargin < 4)
        error('Bad number of arguments');
    end
    %if (nargin < 3)
    %    mask = false(dims.ntotal,1);
    %end
    if (nargin < 3)
        mask = [];
        border = [];
    end
    
    result = [];
    spproj = projector_uniform_unitsimplex(dims);
    result.project = @(u)(project_mask(spproj(u),mask,border)); %SLOW is mask is supplied; many unnecessary operations
    result.center = (1.0 / dims.components) .* ones(dims.ntotal,1);
    % this is the maximal distance in each pixel, i.e.
    % from (1/c,1/c,...,1/c) to (1,0,...,0),
    % times the number of pixels.
    result.radius = sqrt(dims.nimage * ((dims.components-1)/dims.components));
    
    % special functions
    
    result.primal_from_dual = @(v, f, reg)(project_mask(constraints_uniform_unitsimplex_pfd(v, f, reg, dims),mask,border));

    result.support_function = @(x)(supportfunction(x,dims,mask,border));
    
    result.liftable = true;
    result.lift = @(x)([x(:); x(:)]);
    result.unlift = @(x)(sum(reshape(x,[],2),2));    
    result.lifted_project = @(x)(simplex_lift_project(x,dims,mask,border));
    result.lift_factor = 2;

    %result.liftable = true;
    %result.lift = @(x)((1/sqrt(2)) .* [x(:); x(:)]);
    %result.unlift = @(x)((1/sqrt(2)) .* sum(reshape(x,[],2),2));
    %result.lifted_project = @(x)(simplex_lift_project(x,dims));
    %result.lift_factor = 1;

    % Disables lifting
    %result.liftable = true;
    %result.lift = @(x)(x);
    %result.unlift = @(x)(x);
    %result.lifted_project = result.project;
    %result.lift_factor = 1;
    
    % Barrier Function Representation
    
    if (isempty(mask) && isempty(border))
        result.barrier_function.derivatives = @(x)(barrier_derivatives(x,dims));
        result.barrier_function.affine_part = @()(barrier_affine_part(dims));
    end
end

function u = project_mask(u,mask,border)
    u(mask) = border;
end

function result = simplex_lift_project(x, dims, mask, border)
    %SLOW if mask is supplied, many unnecessary operations are done.
    x1 = x(1:dims.ntotal);
    x2 = reshape(x((dims.ntotal+1):end), dims.nimage, dims.components);
    
    % project first vector on x >= 0
    x1 = max(0, x1);
%    s = mean(x2,2) - 1/(sqrt(2) * dims.components);
    s = mean(x2,2) - 1/dims.components;
    x2 = x2 - s(:,ones(1,dims.components));
    
    result = [project_mask(x1(:),mask,border); project_mask(x2(:),mask,border)];
end

function [result,resultlocal] = supportfunction(x,dims,mask,border)
    pfd = project_mask(discretize_firstmax(x(:), dims),mask,border);
    r = pfd(:) .* x(:);
    if (nargout < 2)
        result = sum(r);
    else        
        resultlocal = sum(reshape(r,dims.nimage,dims.components),2);
        result = sum(resultlocal);
    end
%    result = sum(max(reshape(project_mask(x,mask,border), dims.nimage, dims.components),[],2));    
end

% hessianprecond is an approximation of (in fact, in this case exactly)
% (hessian^-1)
function [gradient,hessian,hessianprecond] = barrier_derivatives(x,dims)
    gradient = -1./reshape(x,1,dims.ntotal);
        
    if (nargout > 1)
        hessian = spdiag(gradient.^2);
        
        if (nargout > 2)
            hessianprecond = spdiag(1./(gradient.^2));
        end
    end
end

function [Q,r] = barrier_affine_part(dims)
    [jj,ii] = meshgrid(1:dims.components,1:dims.nimage);

    Q = sparse(...
        flatten(ii), flatten(ii + (jj-1).*dims.nimage), ones(numel(ii),1), ...
        dims.nimage,dims.ntotal); %TODO could also do this via "kron(ones(),speye())"
    r = ones(dims.nimage,1); % note this is non-sparse per default
end

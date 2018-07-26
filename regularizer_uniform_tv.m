%
% Author: Jan Lellmann, j.lellmann@damtp.cam.ac.uk
%

% creates a forward-difference p-norm regularizer with Neumann boundary
% conditions, e.g. lambda * TV(A (*) u), where A (*) u means A is applied
% to u(x) in each point uniformly. If A is not supplied, it is assumed to
% be the identity. A does not necessarily have to be quadratic.
function result = regularizer_uniform_tv(dim, lambda, p, schemes, boundaries, A)   
    if (nargin < 6 || isempty(A))
        A = eye(dim.components); %speye would be nicer but does not work with svd
    end
    components = size(A,1);
    
    righthand = kron(sparse(A),speye(dim.nimage));
    righthandbound = norm(full(A),2); % this should be computable, as A is usually small
    [nabla,nablabound] = diffopn(dim.image, components, schemes, boundaries);

    op = nabla * righthand;
    opnormbound = nablabound * righthandbound; % this is an estimate, but seems to be very tight in practice
   
    % construct dual norm by the rule 1/p+1/dualp=1
    if (p == 1)
        error('p = 1 has no well-defined meaning; consider using regularizer_uniformanisotropic_tv');
        dualp = inf;
        dualradius = sqrt(size(l1,1));
    elseif (p == 2)
        dualp = 2;
        dualradius = sqrt(dim.nimage);
    else
        error('p must be 1 or 2');
    end
    
    result = [];
    
    % primal methods
    result.operator = @(x)(op * x(:));
    result.evaluate = @(x)(regularizer_uniform_tv_evaluate(lambda,reshape(op * x(:), [dim.nimage (dim.mimage * components)]), p));
    result.b = zeros(size(op,1),1);
        
    % dual methods
    result.adjoint_operator = @(y)(op' * y(:));
    
    result.dual_constraints = constraints(...
        @(y)(reshape(project_pnorm(reshape(y,[dim.nimage (dim.mimage * components)]), lambda, dualp),[size(op,1) 1])),... % projector
        lambda * zeros(size(op,1),1),... % center
        lambda * dualradius... % radius %FIXME pessimistic estimate of the RADIUS in 2-norm (here: for dualp = inf)
    );

    result.dual_constraints.liftable = true;
    result.dual_constraints.lift = @(x)(x);
    result.dual_constraints.unlift = @(x)(x);
    result.dual_constraints.lifted_project = result.dual_constraints.project;
    result.dual_constraints.lift_factor = 1; % unlift(lift(x)) = lift_factor * x

    % general information
    result.lengths = size(op);
     
    % set operator norm
    result.operator_norm = opnormbound;        
    
    % special functions    
    result.dual_from_primal = @(u, f)(dfp_regularizer_uniform_tv(u, f, result, lambda, dim));
    
    if (p == 2)
        result.discretization_distance = @(M,N)(lambda *  sqrt(sum(((A*((M-N)'))').^2,2)));
    end
    
    result.operator_backward_step = backward_step_gradient_dct(dim, schemes, boundaries, A); %DEBUG
    if (isempty(result.operator_backward_step))
        result = rmfield(result,'operator_backward_step');
    end 
    
    result.operator_as_matrix = @()(op);
end

function [result,resultlocal] = regularizer_uniform_tv_evaluate(lambda,x,p)
    if (nargout < 2)
        result = lambda * evaluate_pnorm(x,p);
    else
        [result,resultlocal] = evaluate_pnorm(x,p);
        result = result * lambda;
        resultlocal = resultlocal * lambda;
    end        
end

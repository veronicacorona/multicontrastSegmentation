%
% Author: Jan Lellmann, j.lellmann@damtp.cam.ac.uk
%

% implements a variant of the primal-dual method as suggested by Pock et
% al. in "An Algorithm for Minimizing the Mumford-Shah Functional"
% (ICCV09)

% Requires 0 < step_relaxation <= 2
% step_relaxation = 2 is FPD,
% step_relaxation = 1 is exactly the Arrow-Hurwicz/Appleton-Talbot method.

% Notation (here -> original paper)
%   u->x; v->y; ubar->xbar; C->C; D->K; b->h; tau_u->tau; tau_v->sigma;

function [result,details] = solvesfr_primaldual_fpd(start, startdual, dataterm, constraints, reg, params)

    epsilon = getparameter(params,'term_optimality',0); % optimality bound
    relepsilon = getparameter(params,'term_optimality_relative',0); % relative gap (only works well for nonnegative objectives)
    maxiter = ceil(getparameter(params,'term_maxiter',inf)); % optimality bound
    
    callback = getparameter(params,'callback_details',[]); 
    if (isempty(callback))
        callback_plain = getparameter(params,'callback',[]);
        if (isempty(callback_plain))
            callback = [];
        else
            callback = @(u,i,delta,details)(callback_plain(u,i,delta)); % for compatibility
        end
    end

    L = reg.operator_norm;
    tau_u = getparameter(params,'step_primal',0.99 * sqrt(1/(L^2)));
    tau_v = getparameter(params,'step_dual',0.99 * sqrt(1/(L^2))); % Need tau_u * tau_v * L^2 < 1. % TODO could improve weighting
    thetaprime = getparameter(params,'step_relaxation',2.0); % overrelaxation
    
    if (epsilon > 0 || relepsilon > 0)
        if (~dataterm.islinear || ~isfield(constraints,'support_function'))
            error('Evaluation of dual must be supported by constraint set and data term');
        end
        s = dataterm.gradient(start); % get linear part
    end
    
    if (isinf(maxiter) && (epsilon <= 0) && (relepsilon <= 0))
        error('Need at least one valid termination criterion.');
    end   
    
    if (isempty(startdual))
        % choose dual starting point automatically
        startdual = reg.operator(start(:));
    end
    
    % initialization
        
    uk = start(:);
    vk = startdual(:);
    ubark = uk;
    
    if (~isempty(callback))
        details = []; details.dual_feasible_point = vk;
        callback(uk, 0, +inf, details);
    end
    
    for k = 0:(maxiter-1)
        
        % update dual
        
        q = reg.operator(ubark);
        vkp1 = reg.dual_constraints.project(vk + tau_v .* (reg.operator(ubark) - reg.b));
        
        % update primal
        % slightly extended from original paper to use gradient of data term
        % instead of fixed vector
        ukp1 = constraints.project(uk - tau_u .* (reg.adjoint_operator(vkp1) + dataterm.gradient(uk)'));
        
        ubarkp1 = thetaprime * ukp1 + (1-thetaprime) * uk;

        % stop on dual
        if (epsilon > 0 || relepsilon > 0)
            objp = dataterm.evaluate(ukp1) + reg.evaluate(ukp1);
            objd = - constraints.support_function( - reg.adjoint_operator(vkp1(:)) - s(:) );
            gap = objp - objd;
        end

        % ...easy as pie :)
                
        if (~isempty(callback))
            delta = norm(ukp1(:) - uk(:),inf);
            details = [];
            details.dual_feasible_point = vkp1;
            callback(ukp1, k+1, delta, details);
        end
        
        uk = ukp1;
        vk = vkp1;
        ubark = ubarkp1;
        
        if (epsilon > 0)
            if (gap < epsilon)
                break;
            end
        end        
        
        if (relepsilon > 0)
            if ((gap/max(eps,min(abs(objp),abs(objd))) < relepsilon))
                break;
            end
        end        
    end
    
    result = ukp1;
    
    details = [];
    details.dual_feasible_point = vkp1;
    details.iterations = k+1;
end

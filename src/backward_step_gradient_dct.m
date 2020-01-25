%
% Author: Jan Lellmann, j.lellmann@damtp.cam.ac.uk
%

% returns a _function_ @(tau,x) that computes (I + tau*A'*A)^{-1} * x,
% or [] if not possible for given parameters.
function fnhandle = backward_step_gradient_dct(dim, schemes, boundaries, A, scaling)
  if (nargin < 4)
      A = [];
  end
  if (nargin < 5)
      scaling = [];
  end
  fnhandle = [];
  if (strcmp(schemes,'forward') && strcmp(boundaries,'neumann'))
    laplace_factors = laplace_dct_factors_neumann(dim.image, scaling);
    if (isempty(A)) % assume A = I
      %if (size(A,1) == size(A,2) && norm(A - eye(size(A)),'fro') <= eps)
      
      % this stores a copy of laplace_factors in the function pointer's context
      fnhandle = @(tau, x)(reshape(...
          idctmulti(repmat(1./flatten(tau*laplace_factors + 1.0)', dim.components, 1) ...
                         .* dctmulti(flatten(x), dim.image, dim.components),...
                         dim.image, dim.components),[],1));  
%          idctmulti(kron(speye(dim.components),...    
%                        spdiag(1./flatten(tau*laplace_factors + 1.0)))...
%                    * dctmulti(flatten(x), dim.image, dim.components),...
%                    dim.image, dim.components),[],1));
%    end
    else
      k = size(A,1);
      [U,D,R] = svd(A);
      R = R(:,1:k);
      if (isvector(D)) % STUPID Matlab "smart" behavior of diag workaround
          D = D(1,1);
      else
          D = diag(D,0);
      end
      D = D.^2; % the diagonal of A'*A
      
      fnhandle = @(tau, x)(nontrivialaproj(x(:),dim,k,tau,...
          kron(D(:),laplace_factors(:)),R));
    end
  end
end

function result = nontrivialaproj(X,dims,k,tau,DQkron,R)

  %Y = kron(R', speye(dims.nimage)) * X(:); % about 200x slower
  Y = reshape(reshape(X,dims.nimage,dims.components) * R, [], 1);
  
  % DQkron = kron(D(:),Q(:)) is precomputed
  Z = idctmulti((1./(1 + tau.*DQkron)) .* dctmulti(Y,dims.image,k),dims.image,k);
  
  %result = X(:) + kron(R, speye(dims.nimage)) * (-Y + Z); % very slow
  result = X(:) + reshape(reshape(-Y + Z, dims.nimage, k) * R', [], 1);
  
  %Debug code
  %AtA = R*D*R';
  %t = adjop(op(X));  
  %s = kron(R,speye(dims.nimage)) * idctmulti((kron(D,Q)) .* dctmulti(kron(R',speye(dims.nimage))*X,dims.image,k),dims.image,k);
  %norm(s-t,+inf)
  
  %SANITY CHECK
  %disp(['Deviation = ' num2str(norm(result + tau.*adjop(op(result))- X,+inf))]);
end

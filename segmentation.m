% Author: Veronica Corona, vc324@cam.ac.uk

function [fd, u, ud]=segmentation(lambda,k1,post)

dbglevel(5);
rand('seed',192497); % make noise predictable
%load('post3d.mat');
% LOAD LABELS FROM TEMPLATE
nuclei = load_nii('thalamus_all.nii.gz');
n=double(rot90(select_left_thalamus_small(nuclei.img)));

n(n==4)=3;
n(n==5)=1;
nucl=make_labels(n);

if k1==0
    post=post+0.00001;
end

for i = 1:size(post,4)
    f(:,:,:,i)=-log(k1*nucl(:,:,:,i)+(1-k1)*post(:,:,:,i));
end

    dims = dimensions([size(f,1) size(f,2) size(f,3)], 4);


%lambda = 2; % the regularization strength

p = 2; % which norm to use for the regularizer (always use 2 as it is isotropic)
tv = regularizer_uniform_tv(dims, lambda, p, 'forward', 'neumann');
unitsimplex = constraints_uniform_unitsimplex(dims);
dataterm = dataterm_linear(f(:));

info = @(u,i,d,details)print_info(u,i,d,details, 1,dataterm,unitsimplex,tv);
ppd = struct(...
    'term_maxiter',12000,...
    'term_optimality_relative',1e-6,...
    'callback_details',info);

% solve the problem

[u,~] = solvesfr_primaldual_fpd(unitsimplex.center, 0*tv.dual_constraints.center, dataterm, unitsimplex, tv, ppd);

% if we want a discrete choice at each point and no values between two
% labels, we can round u:

ud = discretize_firstmax(u, dims);

% reshape u and ud from vector into image with 4 components (one per
% class)

u = reshape(u, size(f));
ud = reshape(ud, size(f));

fd = discretize_firstmax(-f, dims);
fd = reshape(fd, size(f));

end

function print_info(u,i,~,details, granularity,dataterm,constraints,reg)
if ((mod(i,granularity) == 0))
    obj = dataterm.evaluate(u(:)) + reg.evaluate(u(:));
    objd = - constraints.support_function( - reg.adjoint_operator(details.dual_feasible_point(:)) - dataterm.gradient(constraints.center(:))');
    dbg(1,['#' num2str(i) ', primal = ' num2str(obj) ', dual = ' num2str(objd) ', gap = ' num2str(obj-objd) ', relative gap = ' num2str((obj-objd)/objd)]);
    %        dbg(1,['#' num2str(i)]);
end
end

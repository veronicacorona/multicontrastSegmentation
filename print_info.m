function print_info(u,i,d,details, granularity,dataterm,constraints,reg)
    if ((mod(i,granularity) == 0))
        obj = dataterm.evaluate(u(:)) + reg.evaluate(u(:));
        objd = - constraints.support_function( - reg.adjoint_operator(details.dual_feasible_point(:)) - dataterm.gradient(constraints.center(:))'); 
        dbg(1,['#' num2str(i) ', primal = ' num2str(obj) ', dual = ' num2str(objd) ', gap = ' num2str(obj-objd) ', relative gap = ' num2str((obj-objd)/objd)]);
%        dbg(1,['#' num2str(i)]);
    end
end

tam_prior_eval_log_density_one_parameter <- function( density_pp, args_pp, parameter_pp, eps=1E-20)
{
	args_pp$x <- parameter_pp
	res <- log( do.call( density_pp , args_pp ) + eps )
	return(res)
}
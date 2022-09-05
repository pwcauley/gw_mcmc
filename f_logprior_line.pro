;;General function to calculate a prior log-likelihood for GW_MCMC.PRO based on given input parameters
;; and a given model.  Specifically for F_LINE.PRO

function f_logprior_line,xin,pars

	;;Currently flat prior
	if pars[0] gt -1d5 and pars[0] lt 1d5 and pars[1] gt -10. and pars[1] lt 10. $
		and alog(pars[2]) gt -10. and alog(pars[2]) lt 0. then return,0. $
		else return,-1./0.

end

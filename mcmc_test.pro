;;Test procedure for GW_MCMC.PRO. Will find the best-fit line
;; to a set of points with Gaussian uncertainties, as well as
;;an additional unknown intrinsic scatter that is not accounted
;;for in the known measurement uncertainties. The best-fit is 
;;taken to be the median of each posterior distribution and the
;;1-sigma uncertainties are the standard deviations of the
;;marginalized posteriors.

pro mcmc_test

	;;True line parameters
	b=.86
	y0=-4.2

	;;Generate random points for a line
	npts=50 & sigr=.34 & sig0=.5
	xpts=!RNG->GetRandomNumbers(npts)*10.
	yrand=!RNG->GetRandomNumbers(npts,/normal)
	ypts=y0+xpts*b+yrand*sig0
	yeps=!RNG->GetRandomNumbers(npts,/normal)
	ypts=ypts+ypts*sigr*yeps
	err=fltarr(npts)+sig0
	ltrue=f_line(xpts,[y0,b,sigr])

	;;Set up a plot of the points, the true line, and the best-fit line
	plot,xpts,ypts,ps=2,symsi=1.5,charsi=1.5,/xsty,/ysty,xr=[min(xpts)-1.,max(xpts)+1.],$
		yr=[min(ypts)-1.,max(ypts)+1.],xtit='x-points',ytit='y-points',$
		title='MCMC line fitting test'	
	errplot,xpts,ypts+err,ypts-err

	nwalk0=100 & nsteps0=200
	pg0=[-2.5,1.,.25]
        gw_mcmc,nwalk0,nsteps0,'f_line','f_logprior_line',$
                xpts,ypts,err,pg0,2.,pcout,llout,afout,_EXTRA=extra

	print,'Intercept = '+strtrim(median(pcout(*,0,100:*)),2)+'+-'+strtrim(stddev(pcout(*,0,100:*)),2)
	print,'Slope = '+strtrim(median(pcout(*,1,100:*)),2)+'+-'+strtrim(stddev(pcout(*,1,100:*)),2)
	print,'Intrinsic scatter = '+strtrim(median(pcout(*,2,100:*)),2)+'+-'+strtrim(stddev(pcout(*,2,100:*)),2)

	pmed=[median(pcout(*,0,100:*)),median(pcout(*,1,100:*)),median(pcout(*,2,100:*))]
	lbest=f_line(xpts,pmed)
	oplot,xpts,lbest,co=cgcolor('red'),thick=2
	oplot,xpts,ltrue,co=cgcolor('purple'),thick=2

end

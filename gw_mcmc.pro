;;PROCEDURE GW_MCMC
;;
;;Procedure to perform a Markov Chain Monte Carlo simulation for a specified
;;model and a given data set. The sampling is based on Goodman & Weare's (2010)
;;affine invariant "stretch" method. This program is similar to the Python
;;distribution EMCEE (Foreman-Mackey et al. 2013) except with less general
;;functionality. 
;;
;;Email pwcauley at gmail.com with questions/comments.
;;
;;INPUTS:
;; nwalk - number of walkers
;; nsteps - number of steps in the chain
;; modelfun - the string name of the model function, e.g., 'planetfit'
;; x - vector containing independent variable for data (e.g., wavelength or time)
;; y - vector with same length as X containing the dependent data variable (e.g., flux)
;; err - uncertainty estimates for data Y
;; pguess - initial parameter guesses for the best-fit model
;; afac - scale factor for the density scaling function g(z)
;;        NOTE: This value must be > 1 and the lower it is the
;;              narrower the allowable values of z are from the
;;              scaling function g(z)->[1/z,z]
;; NOTE: Currently ASSUMES that the final parameter is the unknown variance of the data!!
;;	 So PGUESS[npars-1] will be used to estimate the additional data uncertainty.
;;	 This must be a constant fraction of the actual data values. Can still
;;	 be ignored if fitting data where there is no unknown independent variance
;;	 (a high-quality spectrum, for example) 
;;
;;OUTPUTS:
;; pchain - the chain of parameter values with dimension NWALK x NDIM x NSTEPS
;; llchain - the log-likelihood values for each step in the chain for each
;;           walker with dimensions NWALK x NSTEPS
;; afrac - final acceptance fraction of stretch moves
;;
;;KEYWORDS:
;; /loud - prints out progress information
;; ignore= - set equal to an array of locations specifying which
;;           parameters to ignore in the array PGUESS. For example,
;;	     if you want to ignore the second and third parameters in the
;;           parameter vector set ignore=[1,2]
;;           NOTE: The specified parameters are ignored by resetting
;;                 the values to the initial inputs values at each
;;                 step in the chain.
;;
;; parinfo= - string vector containing labels for the parameters, used
;;            for plot labels in MCMC_PLOT.PRO, also useful for ID-ing
;;            parameters upon file restore 
;; /isave - set to save the output, also saves every 100 reps so progress
;;	    isn't lost due to shutdown
;;	    NOTE: If ISAVE is set, you must also specify SDIR=, a string
;;	    variable containing the desired save directory.
;; /timeest - set to print an estimate of the total run-time (not very
;;	      accurate)
;; _EXTRA - Any additional keywords can be specified, by *value* and NOT
;;          by reference, that will be passed to the model function.
;;
;;Written 2016-02-10 by PWC

pro gw_mcmc,nwalk,nsteps,modelfun,priorfun,xin,yin,errin,pguess,afac,$
            pchain,llchain,afrac,loud=loud,ignore=ignore,parinfo=parinfo,$
	    isave=isave,sdir=sdir,timeest=timeest,_EXTRA=extra

	if not keyword_set(isave) then isave=0
	if not keyword_set(ignore) then ignore=-99

	if ignore ne [-99] and ignore eq [-1] then ignore=0

	if n_elements(ignore) gt 0 then begin

		npars=n_elements(pguesS)
		temp=indgen(npars)
		dignore=[]
		for i=0,npars-1 do begin
	
			if where(ignore eq temp[i]) eq [-1] then dignore=[dignore,temp[i]]
	
		endfor
		ndign=n_elements(dignore)

	endif

	;;Prompt user for save directory if not given
	if isave eq 1 and not keyword_set(sdir) then begin

		print,'Save directory SDIR not specified! Please enter string specifying save directory:'
		read,sdir
		
		if strlen(sdir) lt 2 then begin

			print,'SDIR not given in string format, returning...'
			retall

		endif

	endif

	if isave eq 1 then savedir=sdir

  	;;Number of parameters
  	ndim=n_elements(pguess)
	if ignore eq [-99] then ndimq=ndim else ndimq=ndim-n_elements(ignore)

  	if n_elements(parinfo) lt ndim then parinfo=strarr(ndim)+'No info'

  	;;Check for zero errors and don't include these points
  	zz=where(errin eq 0.)
  	if zz ne [-1] then begin
    
		gd=where(errin ne 0.)
    		x=xin(gd) & y=yin(gd) & err=errin(gd)
  
	endif else begin

		x=xin & y=yin & err=errin

	endelse

	;;Save original errors
	err0=err
 
  	;;Initial model guess and log-likelihood
  	if keyword_set(loud) then print,'Calculating initial model...'
	void=execute('mod0='+modelfun+'(x,pguess,_EXTRA=extra)')

	;;Add random uncertainty estimator to errors
	if ignore ne -99 and where(ignore eq ndim-1) eq [-1] then $
		err=sqrt(err0^2.+pguess[ndim-1]^2.*mod0^2.) $
		else err=err0

	log_prior=call_function(priorfun,x,pguess)
  	loglike0=-.5*total(alog(2.*!Pi*err^2.)+((y-mod0)/err)^2.)+log_prior

	;;Flag bad values and restart
	if finite(loglike0) ne 1 then begin

		print,'Initial guesses not allowed by priors! Returning...'
		return

	endif

  	;;Initialize tools for density scaling function g(z). Random deviate z
  	;; chosen based on uniform deviate x and comparing where the area under
  	;; the g(z) curve (approximately) equals the area under the uniform curve f(x).
  	nz=5000
  	z=findgen(nz+1)*afac/float(nz)+1./afac  ;;range of z
  	g=1./sqrt(z)     ;;density function of scaling variable
  	g=g/total(g)
  	areag=[]
  	for i=0,nz-1 do areag=[areag,total(g(0:i))]  ;;area of g as a function of z
  	f=fltarr(nz)+1.
  	f=f/total(f)     ;;uniform distribution over same interval
  
  	;;Initialize walkers. Currently using a normal distribution around
  	;; specified initial parameters with sigma equal to 5% of the 
  	;; parameter value.
  	if keyword_set(loud) then print,'Initializing walkers...'
  	sigmas=.001*pguess
	walk0=!RNG->GetRandomNumbers(nwalk,ndim,/normal)
  	for i=0,ndim-1 do walk0(*,i)=walk0(*,i)*sigmas(i)+pguess(i) 

  	;;IGNORE keyword repeatedly replaces parameter with original value
  	if ignore ne [-99] then begin

    		nign=n_elements(ignore)
    		for ii=0,nign-1 do walk0(*,ignore[ii])=pguess(ignore[ii])

  	endif

  	;;Calculate initial model fits
  	if keyword_set(loud) then print,'Calculating models for initial walker array...'
  	loglike_all=[]
  	for i=0,nwalk-1 do begin

		;;Make sure none of the initial walkers are bad, replace if they are
		repeat begin		

    			log_prior=call_function(priorfun,x,walk0(i,*))

			;;Make sure initial walkers are allowed by priors
			rnums=!RNG->GetRandomNumbers(ndim,/normal)
			if finite(log_prior) ne 1 then walk0(i,*)=pguess+sigmas*rnums

		endrep until finite(log_prior) eq 1 

		void=execute('fit='+modelfun+'(x,walk0(i,*),_EXTRA=extra)')
		if where(ignore eq -99) eq [-1] and where(ignore eq ndim-1) eq -1 then $
			err=sqrt(err0^2.+walk0(i,ndim-1)^2.*fit^2.) $  ;;add in intrinsic scatter
			else err=err0
    		loglike_all=[loglike_all,-.5*total(alog(2.*!Pi*err^2.)+((y-fit)/err)^2.)+log_prior]

  	endfor

  	;;Begin MCMC chain
  	pchain=walk0
  	llchain=loglike_all
  	itot=0. & iacc=0.

        ;;Print estimated time if desired
        if keyword_set(timeest) then tic
 
  	for j=1,nsteps-1 do begin

    		walk_vec=fltarr(nwalk,ndim)
    		loglike_vec=fltarr(nwalk)

    		for i=0,nwalk-1 do begin

      			itot=itot+1.  ;;keep track of total proposed stretch moves

      			walki=pchain(i,*,j-1)  ;;walker that is being tested

      			repeat begin
       
				rnum=!RNG->GetRandomNumbers(1) 
				jj=round(rnum[0]*(float(nwalk)-1.)) ;;choose random X(t) from current chain     
 
			endrep until jj ne i
 
      			cwalk=pchain(jj,*,j-1)   ;;random walker from previous step that is not WALKI
                        rnum=!RNG->GetRandomNumbers(1)
                        xrnd=rnum[0]*(afac-1./afac)+1./afac
      			cprob=total(f(where(z le xrnd)))
      			igz=where(abs(areag-cprob) eq min(abs(areag-cprob)))
      			zfac=z(igz) & zfac=zfac(0)   ;;scale factor for new walker position

      			ynew=cwalk+zfac*(walki-cwalk)    ;;stretch move w/ random walker  

      			;;IGNORE keyword repeatedly replaces parameter with original value
      			if ignore ne [-99] then ynew(ignore)=pguess(ignore)

      			log_prior=call_function(priorfun,x,transpose(ynew))

			if finite(log_prior) eq 1 then begin

				void=execute('modi='+modelfun+'(x,transpose(ynew),_EXTRA=extra)')
				if where(ignore ne -99) eq [-1] and where(ignore eq ndim-1) eq -1 then $
					err=sqrt(err0^2.+ynew[ndim-1]^2.*modi^2.) $
					else err=err0
	      			loglikei=-.5*total(alog(2.*!Pi*err^2.)+((y-modi)/err)^2.)+log_prior

			endif else loglikei=-1./0.

      			qtest=min([0.,float(ndimq-1.)*alog(zfac)+loglikei-llchain(i,j-1)])

			;;Make sure NaNs aren't doing weird things
			if finite(llchain(i,j-1)) ne 1 and finite(loglikei) ne 1 then qtest=-1./0.
                        rtest=alog(!RNG->GetRandomNumbers(1))
                        rtest=rtest[0]

      			if rtest le qtest then begin

        			walk_vec(i,*)=ynew
        			loglike_vec(i)=loglikei
        			iacc=iacc+1.

        			if keyword_set(loud) then begin
 
          				if i eq nwalk-1 then $
              					print,'Acceptance fraction: '+strtrim(float(iacc)/float(itot),2)+$
                    				' at step #'+strtrim(j,2) 

        			endif

      			endif else begin
    
      		  		walk_vec(i,*)=pchain(i,*,j-1)
        			loglike_vec(i)=llchain(i,j-1)  

      			endelse

                        if keyword_set(timeest) and float(j)*float(i) eq 50. and j eq 1 then begin

                                ntest=float(j)*float(i)
                                toc,elapsed_time=etime & etime=etime/ntest
                                print,'Estimated run-time: '+strtrim(etime*float(nwalk)*float(nsteps+1),2)+' seconds'

                        endif

    		endfor    ;;END WALKER LOOP

		pchain=[[[pchain]],[[walk_vec]]]  ;;update chain
    		llchain=[[llchain],[loglike_vec]]

    		;;If doing long chains, intermittently save progress in case of failure
    		if isave eq 1 and j mod 50 eq 0 and j ne 0 then begin

      			afrac=float(iacc)/float(itot)
      			steps_taken=j
			savedir=sdir
      			funname=modelfun

      			save,pchain,llchain,afrac,x,y,err0,parinfo,ignore,steps_taken,funname,$
      				filename=savedir+'mcmc_'+modelfun+'_'+strtrim(ndim,2)+$
               			'_'+strtrim(nsteps,2)+'_'+strtrim(nwalk,2)+'.sav'

    		endif

  	endfor   ;;END STEP LOOP
 
  	afrac=float(iacc)/float(itot)

	funname=modelfun
  	if isave eq 1 then save,pchain,llchain,afrac,x,y,err0,parinfo,ignore,funname,$
      		filename=savedir+'/mcmc_'+modelfun+'_'+strtrim(ndim,2)+$
               '_'+strtrim(nsteps,2)+'_'+strtrim(nwalk,2)+'.sav'

end

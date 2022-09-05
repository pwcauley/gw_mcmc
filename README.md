# gw_mcmc
Markov Chain Monte Carlo procedure similar to the popular Python distribution emcee. Based on the algorithm of Goodman and Weare (2010).
See the header of the procedure for details on calling sequence and inputs/outputs.

This procedure makes use of David Fanning's random number generator object. To make it work you must specify a new system variable, 
preferably in an IDL startup file. You can read about it at David Fannig's website (http://www.idlcoyote.com/code_tips/randomnumbers.html) 
or just put the following line into your IDL startup file:

DefSysV, '!RNG', Obj_New('RandomNumberGenerator')

You can also use the generator by defining the object at the beginning of each session but that's a bit cumbersome, better to have your
startup file do it for you.

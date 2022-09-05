function f_line,xvec,lpars,_EXTRA=extra

	slope=lpars[1]
	yint=lpars[0]

	return,xvec*slope+yint

end

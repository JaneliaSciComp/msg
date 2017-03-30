is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 

if (! is.installed('HiddenMarkov')) {
	cat("R package HiddenMarkov not found\n")
	quit(save='no', status=2)
}

if (! is.installed('R.methodsS3')) {
	cat("R package R.methodsS3 not found\n")
	quit(save='no', status=2)
}

if (! is.installed('R.oo')) {
	cat("R package R.oo not found\n")
	quit(save='no', status=2)
}

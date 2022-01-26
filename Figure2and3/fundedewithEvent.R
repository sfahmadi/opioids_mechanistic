fundedewithEvent<-function(states, fulltimes,truepar,namesyout,eventdata){
	truepar["starttime"]<-unclass(as.POSIXct(strptime(date(),"%c")))[1]
try({out <- dede(states, fulltimes, "derivs", parms=truepar, dllname="delaymymod",initfunc="initmod", nout=33, rtol=1e-14, atol=1e-14, method="adams",n_history = 100000,events=list(data=eventdata))});
colnames(out)[(length(states)+2):(length(states)+length(namesyout)+1)]=namesyout


out
}

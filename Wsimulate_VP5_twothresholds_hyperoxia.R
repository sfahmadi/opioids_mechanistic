require(deSolve)

molar_mass<-336.4

source("fundede.R")
source("fundedewithEvent.R")
source("models/delaystates.R")
source("models/delaypars.R")
#source("models/delaymodelfun.R")
dyn.load("models/delaymymod.so")
set.seed(100)


source("../Events.R"); source("../crossing.R");
#namesyout=c("P_a_co2","C_a_co2","P_a_o2","C_a_o2","C_V_o2","C_V_co2","C_B_o2","C_B_co2","P_A_co2","P_A_o2","C_e_co2","C_e_o2","C_Vb_o2","C_Vb_co2","P_Vb_o2","P_Vb_co2","C_Vt_o2","C_Vt_co2","P_Vt_o2","P_Vt_co2","P_B_co2","Qt","Venti","Plag_f_pc","Plag_P_a_co2","Plag_P_a_co2","Plag_P_a_o2","Qt","Qb","Qb","P_B_o2","Clag_P_B_co2","W-kf")
#imannamesyout=c("P_a_co2","C_a_co2","P_a_o2","C_a_o2","C_V_o2","C_V_co2","C_B_o2","C_B_co2","P_A_co2","P_A_o2","C_e_co2","C_e_o2","C_Vb_o2","C_Vb_co2","P_Vb_o2","P_Vb_co2","C_Vt_o2","C_Vt_co2","P_Vt_o2","P_Vt_co2","P_B_co2","Qt","Venti","Plag_f_pc","Plag_P_a_co2","Plag_P_a_co2","Plag_P_a_o2","psai_co2","psai_o2","Qb","P_B_o2","Clag_P_B_co2","W-kf")
namesyout=c("P_a_co2","C_a_co2","P_a_o2","C_a_o2","C_V_o2","C_V_co2","C_B_o2","C_B_co2","P_A_co2","P_A_o2","C_e_co2","C_e_o2","C_Vb_o2","C_Vb_co2","P_Vb_o2","P_Vb_co2","C_Vt_o2","C_Vt_co2","P_Vt_o2","P_Vt_co2","P_B_co2","Qt","Venti","Plag_f_pc","Plag_P_a_co2","Plag_P_a_co2","Plag_P_a_o2","psai_co2","psai_o2","Qb","P_B_o2","Clag_P_B_co2","W-kf")

#-----------load Experimental data-----------
clinicaldata<-data.frame(time=c(0,20,40,60,120,180,240,300),
		PaCO2=c(39,44,48,50,53,56,59,63),
		PaO2=c(412,423,452,402,385,383,332,314),
		PaO2SD=c(108, 136, 65,16,163,84,93,87))


#imanWbasic=6  
#breath 30 minutes to reach steady state; here ventilation is controlled by ventilator 
Wbasic=5.8      #ventilator setting for baseline ventilation to get ~40 mm Hg baseline PAco2
                #note that ventilation = Wbasic (using the state variable iman) + Dc + alphaH*Dp
                #here because preoxygenation Dp is negative; and for some reason Dc is also negative (CO2 metabolism rate low so brain CO2 level low)
                #so (total) ventilation = Wbasic
Wakefulnessdrive = 0

this.par<-pars
truepar<-this.par[!names(this.par)=="initialdelay" & !names(this.par)=="Dose"] 

#truepar["fL"]=0.04/60/1.2; truepar["fN"]=0.16/60/1.2;

fulltimes<-0:1800
#truepar["offDc"]=0; truepar["offDp"]=0
truepar["W"]= Wakefulnessdrive                               #anesthetization = 0
states["iman"]<-Wbasic                         #ventilator controls basic ventilation
truepar["Bmax"]<-0.66;    #Bmax is the "live space"  
states["P_I_o2"]<-149                        #maintain baseline PaO2 to be >=100 mm Hg
#imanstates["P_I_o2"]<-146
#normal steady state
fulltimes=seq(0,30*60,10)
out=fundede(states=states,fulltimes=fulltimes,truepar=truepar,namesyout=namesyout)
states1=out[nrow(out),names(states)]
#double i_Qb = Qb0*(1+Cim27*(im18+im19))*im25;
#double i_Qt = Qt0*(1+Cim28*rou*im19)*im26;
#then apnea
truepar["Bmax"]=0.66;  #here Bmax is "live space of ventilation"
states1["iman"]<-5.8;    #turn off ventilator
truepar["offO2"]<-1;  #even when airway obstructed there should still be air exchange in lung
truepar["offCo2"]<-1
truepar["Cim27"]<-1
truepar["Cim28"]<-1 #nothing to do in the model 
#truepar["Cim28"]<-1.5
#truepar["Cim28"]<-3 # (Qt+Qb)/(Qt0+Qb0) vs Pao2 is fitted well

fulltimes=seq(0,20*60,0.1)
states1["P_I_o2"]<-760*1/100
#truepar["P_I_co2"]<-0*760

#PRE RUN -----------------------------IMAN------------
Preout1=fundede(states=states1,fulltimes=fulltimes,truepar=truepar,namesyout=namesyout)

adding_threshold_CO2="yes"
adding_threshold_O2="yes"
tsh_value_o2=15
#Iman --differnt criteria finding --------------- 
#tsh_value_co2=53     #No  CA & No recovery 
tsh_value_co2=52     #yes CA & No recovery 
#tsh_value_co2=52.15  #yes CA & yes recovery 
#Iman --differnt criteria finding --------------- 

truepar["tsh_value_o2"]=tsh_value_o2
truepar["tsh_value_co2"]=tsh_value_co2
#using both Co2 and O2 thresholds ---------------------------------
if (adding_threshold_CO2=="yes" & adding_threshold_O2=="yes") {
	P_a_co2_thrsh=tsh_value_co2 #Co2 threshold mmHG
	P_a_o2_thrsh=tsh_value_o2 #O2 thresholds mmHG
}
#using just  O2 thresholds ---------------------------------
if (adding_threshold_CO2=="no" & adding_threshold_O2=="yes") {
	P_a_co2_thrsh=1000 
	P_a_o2_thrsh=tsh_value_o2
}
#using just  CO2 thresholds ---------------------------------
if (adding_threshold_CO2=="yes" & adding_threshold_O2=="no") {
	P_a_co2_thrsh=tsh_value_co2 
	P_a_o2_thrsh=0
}

P_thrsh_O2=Preout1[,c("time","P_a_o2")][Preout1[,c("P_a_o2")]<=P_a_o2_thrsh ,]


P_thrsh_CO2=Preout1[,c("time","P_a_co2")][Preout1[,c("P_a_co2")]>=P_a_co2_thrsh ,]
delaytime=220
if (length(P_thrsh_O2)!=0) { #if we have intersection
truepar["CA_delay_o2"]=P_thrsh_O2[1,"time"]+delaytime
}else{print("Co2 dont cross the threshold")}

if (length(P_thrsh_CO2)!=0) { #if we have intersection
truepar["CA_delay_co2"]=P_thrsh_CO2[1,"time"]+delaytime
}else{print("O2 dont cross the threshold")}


#print(truepar["CA_delay_o2"])
print(c("threshiold is= ",P_a_o2_thrsh,"CA triggerd at=",truepar["CA_delay_o2"]))

#print(truepar["CA_delay_co2"])
print(c("threshiold is= ",P_a_co2_thrsh,"CA triggerd at=",truepar["CA_delay_co2"]))
#----------------------------------------------------
#eventdata_VP5<-data.frame(var=c("D","iman"),time=c(500,500),value=c(0,1),method=c("replace","replace"))
out1=fundede(states=states1,fulltimes=fulltimes,truepar=truepar,namesyout=namesyout)
print("thresholds at real time----------------------")
print(out1[round(out1[,"time"],1)==truepar["CA_delay_o2"],c("P_a_o2")])
print(out1[round(out1[,"time"],1)==truepar["CA_delay_co2"],c("P_a_co2")])

#sss
mergedout<-merge(out1[,c("time","P_a_co2","P_a_o2")],clinicaldata,by="time")
thiserror<- (mergedout[,"P_a_co2"] - mergedout[,"PaCO2"])/mergedout[,"PaCO2"]
errorvec<-thiserror^2

#
#out1=fundede(states=states1,fulltimes=fulltimes,truepar=truepar,namesyout=namesyout)
#mergedout<-merge(out1[,c("time","P_a_co2","P_a_o2")],clinicaldata,by="time")
#thiserror<- (mergedout[,"P_a_co2"] - mergedout[,"PaCO2"])/mergedout[,"PaCO2"]
#errorvec<-thiserror^2

#original model
#truepar["P2"]<-0
#out2=fundede(states=states1,fulltimes=fulltimes,truepar=truepar,namesyout=namesyout)
#colnames(out2)[colnames(out2)=="P_a_co2"]<- "originalP_a_co2"
#colnames(out2)[colnames(out2)=="P_a_o2"]<- "originalP_a_o2"
#mergedout2<-merge(out2[,c("time","originalP_a_co2","originalP_a_o2")],mergedout,by="time")

#plot(out1[,c("time")]/60,out1[,c("Qb")]+out1[,c("Qt")], col="blue", lwd=3,ylab="Cardiac output (l/min)",xlab="time(min)",type="l",cex.lab=1.5,xlim=c(0,20))

out_df<-data.frame(out1)
library(ggplot2)
cols <- c("Clinical"="black","Dynamic"="blue")

p0<-ggplot()
p0<-p0+geom_line(data=out_df,aes(x=time/60,y=(Qb+Qt)*60,col="Dynamic"))
#p0<-p0+geom_point(data=horsedata,aes(x=time,y=PaO2/Baseline_clinO2,color="Clinical"))
p0<-p0+theme_light()
p0<-p0+ theme(legend.position="none")+
		ylab("Cardiac Output (l/min)")+xlab("")
p0<-p0+labs(title=" A   Canine Hypoxia Calibration")
#p0<-p0+geom_vline(xintercept = 4, linetype="dotted", 
#		color = "red", size=1)
#p0<-p0+geom_vline(xintercept = 10, linetype="dotted", 
#		color = "red", size=1)
p0<-p0+scale_colour_manual(name="", values=cols)
p0<-p0+scale_x_continuous(limits=c(0,15))



Canine_data<-read.csv("expdata/Canine_Hypoxia.csv")
p01<-ggplot()
p01<-p01+geom_point(data=Canine_data,aes(x=time,y=MR,color="Clinical"))
p01<-p01+geom_errorbar(data=Canine_data,aes(x=time,ymin=LR,ymax=HR,color="Clinical"))

p01<-p01+theme_light()
p01<-p01+ theme(legend.position="none")+
		ylab("Mean Arterial Pressure (%)")+xlab("")
p01<-p01+labs(title="Canine Clinical Data")
#p01<-p01+geom_vline(xintercept = 4, linetype="dotted", 
#		color = "red", size=1)
#p01<-p01+geom_vline(xintercept = 10, linetype="dotted", 
	#	color = "red", size=1)
p01<-p01+scale_colour_manual(name="", values=cols)
p01<-p01+scale_x_continuous(limits=c(0,15))



pdf("figs/2threshold_iman_VP5_hyperoxia.pdf")
par(mar=c(4,6,4,4)+.1)
#plot(mergedout2[,c("time","PaCO2")],col="red", ylim=c(40,80))
#lines(mergedout2[,c("time","P_a_co2")],col="green")
#lines(mergedout2[,c("time","originalP_a_co2")], col="blue")
#legend("topleft",legend=c("clincial","modified model","original model"),fill=c("red","green", "blue"))
#
##now PaO2
#x<- mergedout2[,"time"]
#y<- mergedout2[,"PaO2"]/mergedout2[1,"PaO2"]
#plot(x,y,col="red", xlab="time (s)", ylab="PaO2 (normalized)", ylim=c(0,1.2))
#epsilon<-10
#sd<-clinicaldata[,"PaO2SD"]
#sd<- sd/clinicaldata[1,"PaO2"]
#segments(x-epsilon,y-sd,x+epsilon,y-sd,col="red")
#segments(x-epsilon,y+sd,x+epsilon,y+sd,col="red")
#segments(x,y-sd,x,y+sd,col="red")
#lines(mergedout2[,"time"], mergedout2[,"P_a_o2"]/mergedout2[1,"P_a_o2"],col="green")
#lines(mergedout2[,"time"], mergedout2[,"originalP_a_o2"]/mergedout2[1,"originalP_a_o2"],col="blue")
#legend("topright",legend=c("clincial"," modified model", "original model"),fill=c("red","green","blue"))
##iman
#
#
#plot(out2[,c("time","originalP_a_co2")], col="blue")
##lines(mergedout2[,c("time","P_a_co2")],col="green")
#lines(out1[,c("time","P_a_co2")], col="green")
#legend("topleft",legend=c("clincial","modified model","original model"),fill=c("red","green", "blue"))

plot(out1[,c("time","P_a_o2")], col="green")
#lines(mergedout2[,c("time","P_a_co2")],col="green")
#lines(out1[,c("time","P_a_co2")], col="green")
#abline(h=.4*out1[1,"P_a_o2"], col="red")
#abline(v=300, col="red")
legend("topleft",legend=c("clincial","modified model","original model"),fill=c("red","green", "blue"))

plot(out1[,c("time","P_a_co2")], col="green")
abline(h=tsh_value_co2, col="red")
#abline(v=300, col="red")
legend("topleft",legend=c("clincial","modified model","original model"),fill=c("red","green", "blue"))


plot(out1[,c("time","P_B_o2")], col="green")
#lines(mergedout2[,c("time","P_a_co2")],col="green")
#lines(out1[,c("time","P_a_co2")], col="green")
legend("topleft",legend=c("clincial","modified model","original model"),fill=c("red","green", "blue"))

plot(out1[,c("time","Qt")], col="green")
lines(out1[,c("time","Qt")], col="green")
abline(h=3.5/4.8*out1[1,"Qt"], col="red")
abline(v=300, col="red")




#lines(mergedout2[,c("time","P_a_co2")],col="green")
#lines(out1[,c("time","P_a_co2")], col="green")
legend("topleft",legend=c("modified model"),fill=c("green"))

plot(out1[,c("time","Qb")], col="green",cex = .5)
lines(out1[,c("time","Qb")], col="green")
abline(v=truepar["CA_delay_o2"], col="red")
abline(v=truepar["CA_delay_co2"], col="blue")
#points(x=truepar["CA_delay_o2"],y=out1[,"Qb"][out1[,"time"]==truepar["CA_delay_o2"]], col="red")
if(truepar["CA_delay_co2"]<5000) {
#points(x=truepar["CA_delay_co2"],y=out1[,"Qb"][out1[,"time"]==truepar["CA_delay_co2"]], col="blue")
}
#lines(out1[,c("time","P_a_co2")], col="green")
legend("topleft",legend=c("modified model"),fill=c("green"))

plot(out1[,c("time")],out1[,c("Qb")]+out1[,c("Qt")], col="green")
abline(h=3.5/4.8*(out1[1,c("Qb")]+out1[1,c("Qt")]), col="red")
abline(v=300, col="red")


plot(out1[,c("time")]/60,out1[,c("Qb")]+out1[,c("Qt")], col="blue", lwd=3,ylab="Cardiac output (l/min)",xlab="time(min)",type="l",cex.lab=1.5,xlim=c(0,20))

#lines(out1[,c("time","Qb")]+out1[,c("time","Qt")], col="green")
#abline(h=3.5/4.8*out1[1,"Qb"], col="red")

#lines(out1[,c("time","P_a_co2")], col="green")
#legend("topleft",legend=c("modified model"),fill=c("green"))


plot(out1[,c("time","yco2")], col="green")
lines(out1[,c("time","psai_co2")],col="red")
#lines(out1[,c("time","P_a_co2")], col="green")
#legend("topleft",legend=c("clincial","modified model","original model"),fill=c("red","green", "blue"))
legend("topleft",legend=c("yco2","psai_co2"),fill=c("green","red"))

plot(out1[,c("time","yo2")], col="green",cex = .5,ylim=c(0,3))
lines(out1[,c("time","psai_o2")],col="red")

#lines(mergedout2[,c("time","P_a_co2")],col="green")
#lines(out1[,c("time","P_a_co2")], col="green")
legend("topleft",legend=c("yo2","psai_o2"),fill=c("green","red"))

plot(out1[,c("time","im25")], col="green")
plot(out1[,c("time","im26")], col="green")

#clinicaldata1<-data.frame(
#		time=c(28.25,32.75,44.75,49.75,79.75),
#		NCardiac=c(1.7,1.4,1.1,1.086538462,1))
#
#
#plot(out1[,"P_a_o2"],(out1[,c("Qb")]+out1[,c("Qt")])/(pars["Qb0"]+pars["Qt0"]), col="green",xlim=c(20,85),ylim=c(0,3),cex = .5)
#points(x=clinicaldata1$time,y=clinicaldata1$NCardiac,xlim=c(20,85))
#
#clinicaldata2<-data.frame(
#
#time=c(21.7721519,26.07594937,25.06329114,26.83544304,31.64556962,33.67088608,39.74683544,53.92405063),
#NCardiac=c(2.616666667,2.5,2.05,2.05,2.016666667,1.75,1.433333333,1.15))
#plot(out1[,"P_a_o2"],(out1[,c("Qb")])/(pars["Qb0"]), col="green",xlim=c(20,85),ylim=c(0,3),cex = .5)
#points(x=clinicaldata2$time,y=clinicaldata2$NCardiac,xlim=c(20,85))


dev.off()
require(deSolve)

molar_mass<-336.4

source("fundede.R")
source("fundedewithEvent.R")
source("modelWithSimpleVentilatoryCollapse/delaystates.R")
source("modelWithSimpleVentilatoryCollapse/delaypars.R")
#source("models/delaymodelfun.R")
dyn.load("modelWithSimpleVentilatoryCollapse/delaymymod.so")
set.seed(100)


source("Events.R"); source("crossing.R");
namesyout=c("P_a_co2","C_a_co2","P_a_o2","C_a_o2","C_V_o2","C_V_co2","C_B_o2","C_B_co2","P_A_co2","P_A_o2","C_e_co2","C_e_o2","C_Vb_o2","C_Vb_co2","P_Vb_o2","P_Vb_co2","C_Vt_o2","C_Vt_co2","P_Vt_o2","P_Vt_co2","P_B_co2","Qt","Venti","Plag_f_pc","Plag_P_a_co2","Plag_P_a_co2","Plag_P_a_o2","psai_co2","psai_o2","Qb","P_B_o2","Clag_P_B_co2","W-kf")

#-----------load Experimental data-----------
clinicaldata<-data.frame(time=c(0,20,40,60,120,180,240,300),
		PaCO2=c(39,44,48,50,53,56,59,63),
		PaO2=c(412,423,452,402,385,383,332,314),
		PaO2SD=c(108, 136, 65,16,163,84,93,87))



#breath 30 minutes to reach steady state; here ventilation is controlled by ventilator 
Wbasic=5.8      #ventilator setting for baseline ventilation to get ~40 mm Hg baseline PAco2
                #note that ventilation = Wbasic (using the state variable iman) + Dc + alphaH*Dp
                #here because preoxygenation Dp is negative; and for some reason Dc is also negative (CO2 metabolism rate low so brain CO2 level low)
                #so (total) ventilation = Wbasic
Wakefulnessdrive = 0

this.par<-pars
truepar<-this.par[!names(this.par)=="initialdelay" & !names(this.par)=="Dose"] 

#truepar["fL"]=0.04/60/1.2; truepar["fN"]=0.16/60/1.2;
#truepar["fL"]=0.04/60/1; truepar["fN"]=0.16/60/1;

fulltimes<-0:1800
#truepar["offDc"]=0; truepar["offDp"]=0
truepar["W"]= Wakefulnessdrive                               #anesthetization
states["iman"]<-Wbasic                         #ventilator controls basic ventilation
truepar["Bmax"]<-0.66;    #Bmax is the "live space"  
states["P_I_o2"]<-600                        #maintain baseline PaO2 to be >=100 mm Hg
#normal steady state
fulltimes=seq(0,30*60,10)
#print(states)
#print("-----------------------------------")
#print(truepar)
out=fundede(states=states,fulltimes=fulltimes,truepar=truepar,namesyout=namesyout)
states1=out[nrow(out),names(states)]
print(head(out))
#then apnea
truepar["Bmax"]=0;  #here Bmax is "live space of ventilation"
states1["iman"]<-0;    #turn off ventilator
truepar["offO2"]<-1;  #even when airway obstructed there should still be air exchange in lung
truepar["offCo2"]<-1
fulltimes=seq(0,35*60,.1)

#print(tail(out))
#print("*****************")
#print(head(out))
#truepar["tau_co2"]<-20
#truepar["tau_o2"]<-10
#--- Add in CPR event

#eventdata<-data.frame(var=c("CPR"),time=c(600),value=c(.25),method=c("replace"))
Preout1=fundede(states=states1,fulltimes=fulltimes,truepar=truepar,namesyout=namesyout)


adding_threshold_CO2="no"
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
}else{print("o2 dont cross the threshold")}

if (length(P_thrsh_CO2)!=0) { #if we have intersection
	truepar["CA_delay_co2"]=P_thrsh_CO2[1,"time"]+delaytime
}else{print("CO2 dont cross the threshold")}


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

#out_1_2["time"]<-seq(700,1200,1)
#out1<-rbind(out,out_1_2)
mergedout<-merge(out1[,c("time","P_a_co2","P_a_o2")],clinicaldata,by="time")
thiserror<- (mergedout[,"P_a_co2"] - mergedout[,"PaCO2"])/mergedout[,"PaCO2"]
errorvec<-thiserror^2
str(out1)
#original model
truepar["P2"]<-0
out2=fundede(states=states1,fulltimes=fulltimes,truepar=truepar,namesyout=namesyout)
colnames(out2)[colnames(out2)=="P_a_co2"]<- "originalP_a_co2"
colnames(out2)[colnames(out2)=="P_a_o2"]<- "originalP_a_o2"
mergedout2<-merge(out2[,c("time","originalP_a_co2","originalP_a_o2")],mergedout,by="time")

O2_canine<-read.csv("expdata/O2_canine.csv",header=TRUE)
CO2_canine<-read.csv("expdata/CO2_canine.csv",header=TRUE)

out_df<-data.frame(out1)
cols <- c("Clinical"="black","Dynamic"="blue")	
library(ggplot2)
p10<-ggplot()
p10<-p10+geom_line(data=out_df,aes(x=time/60,y=P_a_o2,col="Dynamic"))
#p10<-p10+geom_point(data=horsedata,aes(x=time,y=PaO2/Baseline_clinO2,color="Clinical"))
p10<-p10+theme_light()
p10<-p10+ theme(legend.position="none")+
		ylab("Arterial O2 (mmHg)")+xlab("Time (minutes)")
p10<-p10+labs(title="Canine O2 Calibration")
p10<-p10+scale_colour_manual(name="", values=cols)
p10<-p10+scale_y_continuous(limits=c(0,205))


p101<-ggplot()
p101<-p101+geom_point(data=O2_canine,aes(x=time/60,y=MR,col="Clinical"))
#p101<-p101+geom_point(data=horsedata,aes(x=time,y=PaO2/Baseline_clinO2,color="Clinical"))
p101<-p101+theme_light()
p101<-p101+ theme(legend.position="none")+
		ylab("Arterial O2 (mmHg)")+xlab("Time (minutes)")
p101<-p101+labs(title="Canine O2 Clinical")
p101<-p101+scale_colour_manual(name="", values=cols)
p101<-p101+scale_y_continuous(limits=c(0,205))
p101<-p101+scale_x_continuous(limits=c(0,40))



p11<-ggplot()
p11<-p11+geom_line(data=out_df,aes(x=time/60,y=P_a_co2,col="Dynamic"))
#p11<-p11+geom_point(data=horsedata,aes(x=time,y=PaO2/Baseline_clinO2,color="Clinical"))
p11<-p11+theme_light()
p11<-p11+ theme(legend.position="none")+
		ylab("Arterial O2 (mmHg)")+xlab("Time (minutes)")
p11<-p11+labs(title="Canine CO2 Calibration")
p11<-p11+scale_colour_manual(name="", values=cols)
p11<-p11+scale_y_continuous(limits=c(35,105))
p112<-ggplot()
p112<-p112+geom_point(data=CO2_canine,aes(x=time/60,y=MR,col="Clinical"))
#p112<-p112+geom_point(data=horsedata,aes(x=time,y=PaO2/Baseline_clinO2,color="Clinical"))
p112<-p112+theme_light()
p112<-p112+ theme(legend.position="none")+
		ylab("Arterial O2 (mmHg)")+xlab("Time (minutes)")
p112<-p112+labs(title="Canine CO2 Clinical")
p112<-p112+scale_colour_manual(name="", values=cols)
p112<-p112+scale_y_continuous(limits=c(35,105))
p112<-p112+scale_x_continuous(limits=c(0,40))


pdf("figs/VP5_canine_rev01_CA_new.pdf")
par(mar=c(4,6,4,4)+.1)

#plot(out1[,c("time","P_a_co2")], col="green")
plot(out1[,c("time")]/60,out1[,c("P_a_co2")], col="green", lwd=3,ylab="P_a_Co2(mmhg)",xlab="time(min)",type="l",cex.lab=1.5,xlim=c(0,35))

#lines(mergedout2[,c("time","P_a_co2")],col="green")
#lines(out1[,c("time","P_a_co2")], col="green")
#points(CO2_canine[,c("time","MR")],col="black")
#points(CO2_canine[,"time"]/60,CO2_canine[,"MR"],col="black")
points(CO2_canine[,"time"]/60,CO2_canine[,"MR"], col="black", lwd=3,ylab="P_a_Co2 (mmhg)",xlab="time(min)",cex.lab=1.5)

#segments(CO2_canine[,"time"]/60,CO2_canine[,"LR"],CO2_canine[,"time"],CO2_canine[,"HR"],col="black")
#legend("topleft",legend=c("modified model","canine model"),fill=c("green", "black"))

#plot(out1[,c("time","P_a_o2")], col="green")
plot(out1[,c("time")]/60,out1[,c("P_a_o2")], col="green", lwd=3,ylab="P_a_o2(mmhg)",xlab="time(min)",type="l",cex.lab=1.5,xlim=c(0,35))

#lines(mergedout2[,c("time","P_a_co2")],col="green")
#lines(out1[,c("time","P_a_co2")], col="green")
#points(O2_canine[,c("time","MR")],col="black")
#points(O2_canine[,"time"]/60,O2_canine[,"MR"],col="black")
points(O2_canine[,"time"]/60,O2_canine[,"MR"], col="black", lwd=3,ylab="P_a_o2 (mmhg)",xlab="time(min)",cex.lab=1.5)


#segments(O2_canine[,"time"]/60,O2_canine[,"LR"],O2_canine[,"time"],O2_canine[,"HR"],col="black")
#legend("topleft",legend=c("modified model","canine model"),fill=c("green", "black"))

#plot(out1[,c("time","P_B_o2")], col="green")
##lines(mergedout2[,c("time","P_a_co2")],col="green")
##lines(out1[,c("time","P_a_co2")], col="green")
#legend("topleft",legend=c("modified model"),fill=c("green"))
#
#plot(out1[,c("time","Qt")], col="green")
#lines(out1[,c("time","Qt")], col="green")
#segments(850,0,850,.35,col="red")
#segments(0,0,2000,0,"red")
#
#
##lines(mergedout2[,c("time","P_a_co2")],col="green")
##lines(out1[,c("time","P_a_co2")], col="green")
#legend("topleft",legend=c("modified model","LOAF"),fill=c("green","red"))
#
#plot(out1[,c("time","Qb")], col="green")
#lines(out1[,c("time","Qb")], col="green")
##lines(out1[,c("time","P_a_co2")], col="green")
#legend("topleft",legend=c("modified model"),fill=c("green"))
#
#plot(out1[,c("time")],out1[,c("Qb")]+out1[,c("Qt")], col="green")
#abline(h=3.5/4.8*(out1[1,c("Qb")]+out1[1,c("Qt")]), col="red")
#abline(v=300, col="red")
#
##lines(out1[,c("time","Qb")]+out1[,c("time","Qt")], col="green")
##abline(h=3.5/4.8*out1[1,"Qb"], col="red")
#
##lines(out1[,c("time","P_a_co2")], col="green")
#legend("topleft",legend=c("modified model"),fill=c("green"))
#
#
#plot(out1[,c("time","yco2")], col="green")
#lines(out1[,c("time","psai_co2")],col="red")
##lines(out1[,c("time","P_a_co2")], col="green")
##legend("topleft",legend=c("clincial","modified model","original model"),fill=c("red","green", "blue"))
#legend("topleft",legend=c("yco2","psai_co2"),fill=c("green","red"))
#
#plot(out1[,c("time","yo2")], col="green",cex = .5,ylim=c(0,3))
#lines(out1[,c("time","psai_o2")],col="red")
#
##lines(mergedout2[,c("time","P_a_co2")],col="green")
##lines(out1[,c("time","P_a_co2")], col="green")
#legend("topleft",legend=c("yo2","psai_o2"),fill=c("green","red"))
#
#plot(out1[,c("time","im25")], col="green")
#plot(out1[,c("time","im26")], col="green")
dev.off()
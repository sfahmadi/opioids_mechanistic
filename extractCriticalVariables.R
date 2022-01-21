extractCriticalVariables<-function(naloxone_doseN){
	result_Co2TP="result_Co2TP/"
	system(paste0("mkdir -p ",result_Co2TP))		
	Turn_p=c()
	for(timeUL in c(30)){
		for(Nconditions in c(1,2,3,5,6)){#1:5 corresponds to naloxone doseN 0-4 Added 2 doses as per DGS
			ACO2P=pp[[1]][[Nconditions]][,c("time","Arterial CO2 partial pressure (mm Hg)","Blood flow to tissue (l/min)",
							"Blood flow to brain (l/min)","Minute ventilation (l/min)","Arterial O2 partial pressure (mm Hg)",
							"Brain O2 partial pressure (mm Hg)","Brain CO2 partial pressure (mm Hg)")]
			ACO2P[,"time"]<-ACO2P[,"time"]/60 #Convert time to minutes
			n=2.6; k3=26.6; #mm Hg
			PO2Virtual=ACO2P[,"Arterial O2 partial pressure (mm Hg)"]*(40/ACO2P[,"Arterial CO2 partial pressure (mm Hg)"])^0.3;
			SO2=PO2Virtual^n/(k3^n+PO2Virtual^n);
			ACO2P=cbind(ACO2P,oxygenSaturation=SO2*100)
			ACO2P_S=ACO2P[,"oxygenSaturation"]
			idTPS=which(ACO2P_S==min(ACO2P_S))
			idSaturation90=which(abs(ACO2P_S[c(idTPS+1):571]-90)==min(abs(ACO2P_S[c(idTPS+1):nrow(ACO2P)]-90)))			
			dACO2P=diff(ACO2P[,"Arterial CO2 partial pressure (mm Hg)"])
			dAO2P=diff(ACO2P[,"Arterial O2 partial pressure (mm Hg)"])
			dABO2P=diff(ACO2P[,"Brain O2 partial pressure (mm Hg)"])
			dABVP=diff(ACO2P[,"Minute ventilation (l/min)"])
			ACO2P_O=ACO2P[,"Arterial CO2 partial pressure (mm Hg)"]
			ACO2P_B=ACO2P[,"Brain O2 partial pressure (mm Hg)"]
			ACO2P_Ox=ACO2P[,"Arterial O2 partial pressure (mm Hg)"]
			ACO2P_V=ACO2P[,"Minute ventilation (l/min)"]
			idTP=which(diff(sign(dACO2P))!=0, arr.ind=TRUE) #id of turn point
			idTPB=which(ACO2P_B==min(ACO2P_B))
			idTPO2=which(ACO2P_Ox==min(ACO2P_Ox))
			idTPV=which(ACO2P_V==min(ACO2P_V))
			id_Pre_55=which(abs(ACO2P_O[c(idTP+1):nrow(ACO2P)]-55)==min(abs(ACO2P_O[c(idTP+1):nrow(ACO2P)]-55)))
			id_Pre_50=which(abs(ACO2P_O[c(idTP+1):nrow(ACO2P)]-50)==min(abs(ACO2P_O[c(idTP+1):nrow(ACO2P)]-50)))
			id_Pre_45=which(abs(ACO2P_O[c(idTP+1):nrow(ACO2P)]-45)==min(abs(ACO2P_O[c(idTP+1):nrow(ACO2P)]-45)))
			id_Brain_10=which(abs(ACO2P_B[c(idTPB+1):571]-10)==min(abs(ACO2P_B[c(idTPB+1):nrow(ACO2P)]-10)))
			id_Brain_20=which(abs(ACO2P_B[c(idTPB+1):571]-20)==min(abs(ACO2P_B[c(idTPB+1):nrow(ACO2P)]-20)))
			id_Brain_25=which(abs(ACO2P_B[c(idTPB+1):571]-25)==min(abs(ACO2P_B[c(idTPB+1):nrow(ACO2P)]-25)))
			id_O2_25=which(abs(ACO2P_Ox[c(idTPO2+1):571]-25)==min(abs(ACO2P_Ox[c(idTPO2+1):nrow(ACO2P)]-25)))
			id_O2_30=which(abs(ACO2P_Ox[c(idTPO2+1):571]-30)==min(abs(ACO2P_Ox[c(idTPO2+1):nrow(ACO2P)]-30)))
			id_O2_50=which(abs(ACO2P_Ox[c(idTPO2+1):571]-50)==min(abs(ACO2P_Ox[c(idTPO2+1):nrow(ACO2P)]-50)))
			id_O2_65=which(abs(ACO2P_Ox[c(idTPO2+1):571]-65)==min(abs(ACO2P_Ox[c(idTPO2+1):nrow(ACO2P)]-65)))
			id_O2_75=which(abs(ACO2P_Ox[c(idTPO2+1):571]-75)==min(abs(ACO2P_Ox[c(idTPO2+1):nrow(ACO2P)]-75)))
			id_Min_40=which(abs(ACO2P_V[c(idTPV+1):571]-2.7)==min(abs(ACO2P_V[c(idTPV+1):nrow(ACO2P)]-2.7)))
			id_Min_60=which(abs(ACO2P_V[c(idTPV+1):571]-4.1)==min(abs(ACO2P_V[c(idTPV+1):nrow(ACO2P)]-4.1)))
			id_Min_80=which(abs(ACO2P_V[c(idTPV+1):571]-5.4)==min(abs(ACO2P_V[c(idTPV+1):nrow(ACO2P)]-5.4)))
			if (length(idTP)>1) {idTP=idTP[1]}
			Ndosen=(Nconditions-1)
			if ((Nconditions-1)==5) {Ndosen="2+2"}
			Turn_p0=c("opioid"=args$opioid,"patientType"=args$patientType,"dose"=args$conc,"Ndose"=Ndosen,
					ACO2P[idTPB+id_Brain_20,c("time","Minute ventilation (l/min)","oxygenSaturation")],
					ACO2P[idTPO2+id_O2_30,c("time","Minute ventilation (l/min)","oxygenSaturation")],
					ACO2P[idTPB+id_Brain_25,c("time","Minute ventilation (l/min)","oxygenSaturation")],
					ACO2P[idTPO2+id_O2_65,c("time","Minute ventilation (l/min)","oxygenSaturation")],
					ACO2P[idTP+id_Pre_45,c("time","Minute ventilation (l/min)","oxygenSaturation")],
					ACO2P[idTPS+idSaturation90,
							c("time","Minute ventilation (l/min)","Arterial O2 partial pressure (mm Hg)","Arterial CO2 partial pressure (mm Hg)","Brain O2 partial pressure (mm Hg)")]
			)
			Turn_p=rbind(Turn_p,Turn_p0)
		}
	}
	namecsvFile=paste0("opioid=",args$opioid," ","usePK=",args$usePK," ","patientType=",args$patientType," ","dose=",args$conc)
	write.csv(Turn_p,paste0(result_Co2TP,namecsvFile,".csv"),row.names=F)
}
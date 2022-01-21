source("Wsimulate_VP5_canine.R")
print("AAAAs")
source("Wsimulate_VP5_twothresholds_apnea_horse.R")
print("BBBBs")
source("Wsimulate_VP5_twothresholds_hyperoxia.R")

print("CCCCs")
library(ggplot2)
library(grid)
library(gridExtra)
p1<-ggplot()

#pall1<-grid.arrange(p4,p44,p5,p55,
#		p6,p1, p10,p101,
#		p11,p112,
#		ncol=2)
blank <- grid.rect(gp=gpar(col="white"))

pall1<-grid.arrange(p4,p44,p5,p55,ncol=2)
pall2<-grid.arrange(p6,blank,ncol=2)
pall3<-grid.arrange(p10,p101,p11,p112,ncol=2)
pall4<-grid.arrange(p7,p8,ncol=2)
pall5<-grid.arrange(p6,blank,p0,blank,ncol=2)

p0<-p0+labs(title="A Calibration of Cardiac Output")

p01<-p01+labs(title=" Canine Cardiac Output",y="Mean Arterial Pressure (%control)")
p6<-p6+labs(title="B Calibration of Brain Blood Flow")
p66<-p66+labs(title = "Porcine Brain Blood Flow", y = "CBF (ml/100g/min)")
p4<-p4+labs(title="C PaO2 Validation", y="PaO2 (Fraction of control)")
p5<-p5+labs(title="PaCO2 Validation",y="PaCO2 (Fraction of control)")
p44<-p44+labs(title="Equine PaO2")
p55<-p55+labs(title="Equine PaCO2")
pall6<-grid.arrange(p0,p01,p6,p66,p4,p44,p5,p55,ncol=2)
ggsave("Manuscript_figs/Calibration_Validation_Digitized3.pdf",pall6,width=8,height=10.5)

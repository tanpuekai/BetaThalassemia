library(ggplot2)
library(dplyr)
library(grid)
library(ggbeeswarm)
library(png)
library(gridExtra)

BTClinical<-read.table("clinical-311-BT-patients-apr15-2021.txt",header=T,sep="\t")
BTClinical.Apr<-readr::read_tsv("clinical-311-BT-patients-apr15-2021.txt")
BTClinical<-readr::read_tsv("clinical-311-BT-patients-May6-2021.txt")
RNAseq<-readr::read_tsv("clinical-for-RNAseq.with.names.May6.txt")
RNAseq2<-RNAseq%>%group_by(NameChi)%>%slice_head(n=1)%>%arrange(PatientIDs) %>% 
	mutate(Ferro=as.numeric(SerFer),log10Ferro=log10(Ferro),EngName=NameChi,FerGrading=BTGroupChi)

setdiff(RNAseq2$NameChi,BTClinical$ChiName)
ClinRNA<-BTClinical[match(RNAseq2$NameChi,BTClinical$ChiName),]

table(ClinRNA$FerGrading0,RNAseq2$BTGroup)
sapply(split(RNAseq2$NameChi,paste0(ClinRNA$FerGrading0,";",RNAseq2$BTGroup)),paste0,collapse=",")

table(BTClinical.Apr$ChiName==BTClinical$ChiName)
table(BTClinical.Apr$Gender==BTClinical$Gender)
table(BTClinical.Apr$ALP==BTClinical$ALP)
table(BTClinical.Apr$AgeStr==BTClinical$AgeStr)
################################
Yr<-as.integer(gsub("year.*","",BTClinical.Apr$AgeStr))
Mo<-as.integer(gsub("month","",gsub(".*year|[0-9]+day$","",BTClinical.Apr$AgeStr)))
Dy<-as.integer(gsub("day","",gsub(".*month","",gsub(".*year","",BTClinical.Apr$AgeStr))))
Mo[is.na(Mo)]<-0;Dy[is.na(Dy)]<-0;

Age<-Yr+Mo/12+Dy/365
BTClinical.Apr[["Age"]]<-Age

plot(BTClinical.Apr$Age,BTClinical$Age)

table(BTClinical.Apr$preTranspFer==BTClinical$preTranspFer)
plot(as.numeric(BTClinical.Apr$preTranspFer),as.numeric(BTClinical$preTranspFer))
data.frame(BTClinical.Apr$preTranspFer[which(BTClinical.Apr$preTranspFer!=BTClinical$preTranspFer)],
	BTClinical$preTranspFer[which(BTClinical.Apr$preTranspFer!=BTClinical$preTranspFer)])

table(BTClinical$Ferro<3e3 & BTClinical$Age<4,BTClinical$FerGrading)
table(BTClinical$Ferro>5e3 & BTClinical$Age>8,BTClinical$FerGrading)

################################
Fer<-as.numeric(gsub(">","",BTClinical$preTranspFer))
hist(Fer,breaks=1e2)
qqnorm(log(Fer))

BTClinical[["Ferro"]]<-Fer

BTClinical[["log10Ferro"]]<-log10(Fer)


fitAge<-lm(log10Ferro~Age,BTClinical)
boxplot(fitAge$residuals~BTClinical$Gender[!is.na(BTClinical$log10Ferro)])
summary(lm(fitAge$residuals~BTClinical$Gender[!is.na(BTClinical$log10Ferro)]))

str1<-capture.output(summary(lm(log10Ferro~Gender+Age,BTClinical)))
str2<-capture.output(summary(lm(log10Ferro~Gender+Age+ALP,BTClinical)))

str3<-capture.output(summary(lm(ALP~Gender+Age,BTClinical)))
str4<-capture.output(summary(lm(ALP~Gender+Age,BTClinical)))

BTClinical[["FerroResidualsFittedOnAge"]]<-fitAge$residuals[match(seq_along(BTClinical$index),names(fitAge$residuals))]
################################
mytheme2 <- gridExtra::ttheme_default(padding = unit(c(1, 1), "mm"),
    		core = list(fg_params=list(cex = 2/5)),
    		colhead = list(fg_params=list(cex =  2/5)),
    		rowhead = list(fg_params=list(cex =  2/5)))
################################

pdf("BT-ferro-ALP-gender-analyses.311.patients.May6.pdf",width=10)
	grid.arrange(textGrob(paste0("This is a new analyses on ", date()," by Dr Peikai Chen
on the Iron Overload (IO) issue in 311 BT-major patients, where
the I,II,III grading are now provided.

We analysed the ferritin vs age, liver size, or ALP.")))

	myt0 <- gridExtra::tableGrob(BTClinical%>%slice_head(n=70), theme=mytheme2 )
	grid.arrange(myt0,top="Summary of the top 70 patients (total 311 patients) in current analyses")

####

	grid.arrange(textGrob("First, let's take a look at data distributions."))
	p1<-ggplot(BTClinical,aes(x=Ferro)) +
		geom_histogram(aes(y=..density..),binwidth=250) +
		geom_density(alpha=.2, fill="#FF6666") +
		ggtitle("Histogram of Ferro distributions among the 311 patients (308 have values, 3 don't)")

	p2<-ggplot(BTClinical,aes(x=log10(Ferro))) +
		geom_histogram(aes(y=..density..)) +
		geom_density(alpha=.2, fill="#FF6666") +
		ggtitle("Histogram of Ferro distributions among the 311 patients (308 have values, 3 don't)")
	grid.arrange(p1,p2)

	qqnorm(scale(log10(Fer)))
	abline(0,1,col=2)
	legend("topleft",legend="This chart shows that the Ferro is perfectly Gaussian after taking logarithm, 
which means for whatever subsequent analyses, 
we need to be based on the logarithmic scale.",bty='n')

####
	ObChi<-ggplot(BTClinical,aes(x=Age,y=log10(Ferro),label=ChiName,group=Gender)) +
		geom_smooth(data=BTClinical%>%filter(Gender=="Male"),aes(x=Age,y=log10(Ferro)),color="blue") +
		geom_smooth(data=BTClinical%>%filter(Gender=="Female"),aes(x=Age,y=log10(Ferro)),color="red") +
		geom_point(aes(color=Gender),size=2) +
		geom_text() +
		ggtitle("Serum ferritin concentration vs. age and gender in 311 BT patients.
Regression based on LOESS method.") +
		ylab(expression(log[10](Serum~Ferritin~Conc.))) + xlab("Age (Yrs)")
	ggsave("ChiName.png",width=10)

	pngRaster <- readPNG("ChiName.png")
	grid.newpage()
 	grid.raster(pngRaster, width=unit(0.98, "npc"), height= unit(0.98, "npc"))

####
	ObChi<-ggplot(BTClinical,aes(x=Age,y=log10(Ferro),label=ChiName,group=Gender)) +
		geom_smooth(data=BTClinical%>%filter(Gender=="Male"),aes(x=Age,y=log10(Ferro)),color="blue",method="lm") +
		geom_smooth(data=BTClinical%>%filter(Gender=="Female"),aes(x=Age,y=log10(Ferro)),color="red",method="lm") +
		geom_point(aes(color=Gender),size=2) +
		geom_text() +
		ggtitle("Serum ferritin concentration vs. age and gender in 311 BT patients
	(Grouped by gender)
Regression based on linear method.") +
		ylab(expression(log[10](Serum~Ferritin~Conc.))) + xlab("Age (Yrs)")
	ggsave("ChiName.png",width=10)

	pngRaster <- readPNG("ChiName.png")
	grid.newpage()
 	grid.raster(pngRaster, width=unit(0.98, "npc"), height= unit(0.98, "npc"))


####

	ObChi<-ggplot(BTClinical,aes(x=Age,y=log10(Ferro),label=EngName,group=FerGrading)) +
		geom_smooth(data=BTClinical%>%filter(FerGrading0=="￠?"),aes(x=Age,y=log10(Ferro)),method="lm",color="blue") +
		geom_smooth(data=BTClinical%>%filter(FerGrading0=="￠ò"),aes(x=Age,y=log10(Ferro)),method="lm",color="red") +
		geom_smooth(data=BTClinical%>%filter(FerGrading0=="￠ó"),aes(x=Age,y=log10(Ferro)),method="lm",color="green") +
		geom_point(aes(color=FerGrading0,shape=FerGrading0),size=4) +
		ggtitle("Serum ferritin concentration vs. age
in a gender-specific manner, in 311 BT patients
	Using Linear Regression
	Color-coded by FER-grading (your codes)") +
		ylab(expression(log[10](Serum~Ferritin~Conc.))) + xlab("Age (Yrs)")


	ggsave("ChiName.png",width=10)
	pngRaster <- readPNG("ChiName.png")
	grid.newpage()
 	grid.raster(pngRaster, width=unit(0.98, "npc"), height= unit(0.98, "npc"))

####
	ObChi<-ggplot(BTClinical,aes(x=Age,y=log10(Ferro),label=EngName,group=FerGrading)) +
		geom_smooth(data=BTClinical%>%filter(FerGrading0=="￠?"),aes(x=Age,y=log10(Ferro)),method="lm",color="blue") +
		geom_smooth(data=BTClinical%>%filter(FerGrading0=="￠ò"),aes(x=Age,y=log10(Ferro)),method="lm",color="red") +
		geom_smooth(data=BTClinical%>%filter(FerGrading0=="￠ó"),aes(x=Age,y=log10(Ferro)),method="lm",color="green") +
		geom_point(aes(color=FerGrading0,shape=FerGrading0),size=4) +
		geom_text(data=BTClinical%>%filter(FerGrading0%in%c("￠?","￠ó")),aes(x=Age,y=log10(Ferro),label=ChiName)) +
		ggtitle("Serum ferritin concentration vs. age 
		in a gender-specific manner, in 311 BT patients
	Using Linear Regression
	Color-coded by FER-grading (your codes) and labelled with Chinese names") +
		ylab(expression(log[10](Serum~Ferritin~Conc.))) + xlab("Age (Yrs)")
	ggsave("ChiName.png",width=10)
	pngRaster <- readPNG("ChiName.png")
	grid.newpage()
 	grid.raster(pngRaster, width=unit(0.98, "npc"), height= unit(0.98, "npc"))

###############

	ObChi<-ggplot(BTClinical,aes(x=Age,y=log10Ferro,label=EngName,group=FerGrading)) +
		geom_point(aes(color=FerGrading0),size=4) +
		geom_point(data=RNAseq2,mapping=aes(as.numeric(Age),log10Ferro,color=FerGrading),size=6)+
		geom_text(data=RNAseq2,mapping=aes(as.numeric(Age),log10Ferro),size=6)+
		ggtitle(paste0("Serum ferritin concentration vs. age 
		in a gender-specific manner, in 311 BT patients, plus ",paste0(setdiff(RNAseq2$NameChi,BTClinical$ChiName)[1:2],collapse=" and "),"
	Using Linear Regression
	Color-coded by FER-grading (your codes).
Sample labelled with Chinese names are RNAseq samples.")) +
		ylab(expression(log[10](Serum~Ferritin~Conc.))) + xlab("Age (Yrs)")

	ggsave("ChiName.png",width=10)
	pngRaster <- readPNG("ChiName.png")
	grid.newpage()
 	grid.raster(pngRaster, width=unit(0.98, "npc"), height= unit(0.98, "npc"))
###############
	tab12<-sapply(unique(paste0(ClinRNA$FerGrading0)),function(x){
		sapply(unique(paste0(RNAseq2$FerGrading)),function(y){
			ov1<-intersect(ClinRNA$ChiName[ClinRNA$FerGrading0==x],
				RNAseq2$NameChi[RNAseq2$FerGrading==y])
			if(length(ov1)>0)
				paste0(length(ov1),"\n(",paste0(ov1,collapse=", "),")")
			else
				return(0)
		})
	})
	png("Yaya.vs.DrQu.png",1200,1280)
	grid.draw(tableGrob(tab12[c(2,3,4,1),c(1,4,3,2)]))
	grid.draw(textGrob("Summary of the grading of the RNAseq samples, by Yaya and Dr Qu\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"))
	dev.off()

	pngRaster <- readPNG("Yaya.vs.DrQu.png")
	grid.newpage()
 	grid.raster(pngRaster, width=unit(0.98, "npc"), height= unit(0.98, "npc"))
###############
	myt0 <- gridExtra::tableGrob(RNAseq2)
	myt1 <- gridExtra::tableGrob(ClinRNA)
	png("RNAseq.samples.png",1200,1280)
	grid.arrange(myt0,myt1 , top="Summary of the RNAseq samples")
	dev.off()

	pngRaster <- readPNG("RNAseq.samples.png")
	grid.newpage()
 	grid.raster(pngRaster, width=unit(0.98, "npc"), height= unit(0.98, "npc"))


###############
	grid.arrange(textGrob("Next, we investigate Liver size with respect to Age and Ferro.

	(We only have 71 records with Liver size values, of which only
37 are greater than zero. See next page.)"))

myt0 <- gridExtra::tableGrob(BTClinical%>%filter(!is.na(Liver )), theme=mytheme2 )
	grid.arrange(myt0,top="Summary of the 71 patients with Liver size info")


###############
	ggplot(BTClinical,aes(FerGrading0,Liver))+
		geom_violin()+geom_quasirandom(aes(color=FerGrading0),size=4)+
		ggtitle("Liver size vs age, color-coded by grades")

###############

	p1<-ggplot(BTClinical,aes(Liver,log10Ferro))+
		geom_smooth(method="lm",color="blue") +
		geom_point(aes(color=FerGrading0),size=4)+
		ggtitle(paste0("Liver size vs log10Ferro, color-coded by grades\n\n",
			paste0(capture.output(summary(lm(log10Ferro~Liver,BTClinical)))[-seq(7)],collapse="\n")))+
		   	theme(plot.title = element_text(size = 8, family="mono"))

	p2<-ggplot(BTClinical,aes(Liver,log10Ferro))+
		geom_smooth(method="lm",color="blue") +
		geom_quasirandom(aes(color=FerGrading0),size=4)+
		ggtitle(paste0("Liver size vs log10Ferro, color-coded by grades\n(same as left plot, but dots jittered)",
			paste0(rep("\n",12),collapse="")))+
		   		theme(plot.title = element_text(size = 8))
	grid.arrange(p1,p2,ncol=2)
###############
	ggplot(BTClinical,aes(Age,Liver))+
		geom_smooth(method="lm",color="blue") +
		geom_point(aes(color=FerGrading0),size=4)+
		ggtitle(paste0("Liver size vs age, color-coded by grades\n\n",
			paste0(capture.output(summary(lm(Liver~Age,BTClinical)))[-seq(7)],collapse="\n")))+
		   	theme(plot.title = element_text(size = 8, family="mono"))

###############
	t0<-textGrob("Liver size does not depend on Ferro. It depends only on age.")
	t1<-textGrob(paste0(capture.output(summary(lm(log10Ferro~Liver+Age,BTClinical))),collapse="\n"),  gp = gpar(fontsize=9,fontfamily="mono"))
	t2<-textGrob(paste0(capture.output(summary(lm(Liver~log10Ferro+Age,BTClinical))),collapse="\n"),  gp = gpar(fontsize=9,fontfamily="mono"))

	grid.arrange(t0,t1,t2,layout_matrix=matrix(c(1,1,2,3),byrow=T,ncol=2),heights=c(1,8))
	grid.lines(c(0.5,0.5),c(0,1))
	grid.lines(c(0,1),c(0.8,0.8))

	DatLiv<-data.frame(BTClinical%>%filter(!is.na(BTClinical$Liver)),AgeCorrectedLiverSize=lm(Liver~Age,BTClinical)$residuals)

###############

	ggplot(DatLiv,
		aes(log10Ferro,AgeCorrectedLiverSize)) + geom_smooth(method="lm") +
		geom_point(aes(color=FerGrading),size=4)+xlab(expression(log[10](Serum~Ferritin)))+
		ggtitle(paste0("Liver size vs logFerro, color-coded by grades\n\n",
			paste0(capture.output(summary(lm(AgeCorrectedLiverSize~log10Ferro,DatLiv))),collapse="\n")))+
		   	theme(plot.title = element_text(size = 8, family="mono"))

	grid.arrange(textGrob("Basically, we can conclude Liver size does not depend on Ferro, 
nor the other way round. It only depends on age. 

So now we know both Liver size and ALP are irrelevant to serum ferro.

Also Ferro heavily depends on age. So next we corrected the Ferro by the age,
and compare it against gender, grading (I, II, III)"))

###############
	ggplot(BTClinical,aes(x=Age,y=FerroResidualsFittedOnAge,label=EngName,group=Gender)) +
		geom_smooth(data=BTClinical%>%filter(Gender=="Male"),aes(x=Age,y=FerroResidualsFittedOnAge),method="lm",color="blue") +
		geom_smooth(data=BTClinical%>%filter(Gender=="Female"),aes(x=Age,y=FerroResidualsFittedOnAge),method="lm",color="red") +
		geom_point(aes(color=Gender),size=4) +
		ggtitle("Serum ferritin concentration vs. age and gender in 311 BT patients
	Using Linear Regression \n(using the residuals corrected for age)") +
		ylab(expression(log[10](Serum~Ferritin~Conc.~Age~Residuals))) + xlab("Age (Yrs)")


	ggplot(BTClinical,aes(x=Age,y=FerroResidualsFittedOnAge,label=EngName,group=Gender)) +
		geom_smooth(data=BTClinical%>%filter(Gender=="Male"),aes(x=Age,y=FerroResidualsFittedOnAge),method="lm",color="blue") +
		geom_smooth(data=BTClinical%>%filter(Gender=="Female"),aes(x=Age,y=FerroResidualsFittedOnAge),method="lm",color="red") +
		geom_point(aes(color=FerGrading),size=4) +
		ggtitle("Serum ferritin concentration vs. age and gender in 311 BT patients
	Using Linear Regression \n(using the residuals corrected for age), colorcoded by BT grading") +
		ylab(expression(log[10](Serum~Ferritin~Conc.~Age~Residuals))) + xlab("Age (Yrs)")

	grid.arrange(tableGrob(round(100*sapply(split(rank(BTClinical$FerroResidualsFittedOnAge),BTClinical$FerGrading),summary)/nrow(BTClinical),2)),
		top="Summary of percentiles for the three grades \n(BT-1/2/3, corresponding to BT-I/II/III, respectively)")

###############
	ggplot(BTClinical,aes(x=FerGrading,y=FerroResidualsFittedOnAge)) +
		geom_violin(aes(color=FerGrading),size=2) + geom_boxplot(width=0.5) + geom_quasirandom(aes(color=FerGrading),size=4) +
		ggtitle("Age corrected serum ferritin concentration vs. BT grades in 311 BT patients
	(BT-1/2/3, correspond to BT-I/II/III, respectively)")


###############
	str1<-capture.output(summary(lm(log10Ferro~Gender+Age,BTClinical)))
	fitAge<-lm(log10Ferro~Age,BTClinical)
	str2<-capture.output(summary(lm(FerroResidualsFittedOnAge~Gender+Age,BTClinical)))
	str1[3]<-gsub(".*mula = |\\)","",str1[3])
	str2[3]<-gsub(".*mula = |\\)","",str2[3])
	tG1 <- textGrob(paste0(str1[-c(1:2)],collapse="\n"), just="center",  gp = gpar(fontsize=9,fontfamily="mono"))
	tG2 <- textGrob(paste0(str2[-c(1:2)],collapse="\n"), just="center",  gp = gpar(fontsize=9,fontfamily="mono"))

	grid.arrange(tG1,tG2)
	grid.lines(c(0,1),c(0.5,0.5))

	str0<-capture.output(summary(lm(log10Ferro~Age,BTClinical)))
	str0[3]<-gsub(".*mula = |\\)","",str0[3])
	tG0 <- textGrob(paste0(str0[-c(1:2)],collapse="\n"), just="center",  gp = gpar(fontsize=9,fontfamily="mono"))
	tG03<-textGrob(expression(atop("What this regression means is that:"~ Ferro == 10^{3.32+0.025*italic(x)} ,
		 "or " ~ Ferro == 2089.296 %*% 1.059254^italic(x)~"where" ~italic(x)~ "is age in years")))
	grid.arrange(tG0,tG03)
	grid.lines(c(0,1),c(0.5,0.5))


	ggplot(BTClinical,aes(x=Gender,y=FerroResidualsFittedOnAge,label=EngName,group=Gender)) +
		geom_violin() +
		geom_boxplot(width=0.5) +
		geom_quasirandom(aes(color=Gender),size=2) +
		ggtitle("After corrected for age, the residuals show no \nsignificant difference betwn the genders") +
		ylab(expression(Residuals~of~log[10](Ferro)~Fitted~On~Age)) + xlab("Gender") 

	MatVAR<-cbind("std dev"=signif(sd(BTClinical$FerroResidualsFittedOnAge,na.rm=T)),
				"variance"=signif(var(BTClinical$FerroResidualsFittedOnAge,na.rm=T)))

	ggplot(BTClinical,aes(x=FerroResidualsFittedOnAge)) +
		geom_histogram(aes(y=..density..)) +
		geom_density(alpha=.2, fill="#FF6666") +
		ggtitle("Histogram of Corrected Ferro distributions \namong the 311 patients (308 have values, 6 don't)

We define the lower 25% as mild, middle 50% as moderate, and upper 25% as severe") +
		xlab("Corrected Ferro (corrected by age)")+
		annotation_custom(tableGrob(MatVAR), 
          	xmin = -0.8, xmax = -0.4, ymin = 1.0, ymax = 2.0) +
		geom_vline(xintercept=quantile(BTClinical$FerroResidualsFittedOnAge,c(1/3,2/3),na.rm=T))

	fit2<-lm(log10Ferro~Age,BTClinical)

	ageadjust<-coefficients(fit2)[2]*seq(20)+coefficients(fit2)[1]
	ageadjust_min<-ageadjust + quantile(BTClinical$FerroResidualsFittedOnAge,1/3,na.rm=T)
	ageadjust_max<-ageadjust + quantile(BTClinical$FerroResidualsFittedOnAge,2/3,na.rm=T)

	DATREF<-data.frame(years=seq(20), ferro_min=round(10^(ageadjust_min),2),
		ferro=round(10^(ageadjust),2),
		ferro_max=round(10^(ageadjust_max),2))

	grid.arrange(tableGrob(DATREF))
#######################################
	grid.arrange(textGrob("So according to this new metric, the RNAseq samples would look like this (next page)"))

	IndGrading<-sapply(seq_along(RNAseq2$Age),function(i){
		thisAge<-as.numeric(RNAseq2$Age)[i]
		
		indAGE<-which.min(abs(thisAge-seq(20)))
		return(indAGE[1])
	})
	newgrading<-rep("Medium",nrow(RNAseq2))
	newgrading[which(RNAseq2$Ferro<DATREF$ferro_min[IndGrading])]<-"Low"
	newgrading[which(RNAseq2$Ferro>DATREF$ferro_max[IndGrading])]<-"High"
	newgrading[grep("CK",RNAseq2$PatientIDs)]<-"Controls"

	RNAseq2[["newgrading"]]<-newgrading
	ClinRNA[["newgrading"]]<-newgrading



	myt0 <- gridExtra::tableGrob(RNAseq2%>%select(!c(Gender,Age,SerFer)))
	myt1 <- gridExtra::tableGrob(ClinRNA%>%select(!c(ALP,ageYear, ageMonth ,ageDay, EngName,preTranspFerChi)))
	png("RNAseq.samples.png",1200,1280)
	grid.arrange(myt0,myt1 , top="If we use our new grading ...")
	dev.off()

	pngRaster <- readPNG("RNAseq.samples.png")
	grid.newpage()
 	grid.raster(pngRaster, width=unit(0.98, "npc"), height= unit(0.98, "npc"))
####

	png("RNAseq.samples.new.greading.png",1200,1280)
	grid.arrange(tableGrob(table(newgrading,RNAseq2$FerGrading)[c(1,3,4,2),c(1,3,2,4)]),
		tableGrob(table(newgrading,ClinRNA$FerGrading)[c(1,3,4,2),]),
		top="and these are the two confusion tables")
	dev.off()

	pngRaster <- readPNG("RNAseq.samples.new.greading.png")
	grid.newpage()
 	grid.raster(pngRaster, width=unit(0.98, "npc"), height= unit(0.98, "npc"))
###############
	tab12<-t(sapply(unique(paste0(RNAseq2$newgrading)),function(x){
		sapply(unique(paste0(RNAseq2$FerGrading)),function(y){
			ov1<-intersect(RNAseq2$NameChi[RNAseq2$newgrading==x],
				RNAseq2$NameChi[RNAseq2$FerGrading==y])
			if(length(ov1)>0)
				paste0(length(ov1),"\n(",paste0(ov1,collapse=", "),")")
			else
				return(0)
		})
	}))[c(3,4,2,1),c(2,3,4,1)]

	png("Yaya.vs.Peikai.png",1200,1280)
	grid.draw(tableGrob(tab12))
	grid.draw(textGrob("Summary of the grading of the RNAseq samples, by Yaya and Peikai\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"))
	dev.off()

	pngRaster <- readPNG("Yaya.vs.Peikai.png")
	grid.newpage()
 	grid.raster(pngRaster, width=unit(0.98, "npc"), height= unit(0.98, "npc"))
dev.off()

save(RNAseq2,file="NewGrading.Peikai.RData")




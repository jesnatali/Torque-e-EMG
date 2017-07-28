rm(list=ls())
# Download dos pacotes exigidos - Necessario instalar apenas uma vez por R
install.packages("ggplot2")
install.packages("reshape2")
install.packages("readxl")
#install.packages("biosignalEMG")
install.packages("signal")
install.packages("caTools")
install.packages("gridExtra")


# Habilitar pacotes exigidos
#   Caso nao estejam instalados, executar as linhas comentadas acima
library(ggplot2)
library(reshape2)
library(readxl)
library(tools)
#library(biosignalEMG)
library(signal)
library(caTools)
library(gridExtra)


# Endereco das funcoes R para executar o procedimento do ECG Relatorio - Windows
source(file="/Users/Du/Dropbox/PosDoc/Projetos paralelos/EMG e Torque 2/FuncoesForca.R")

# Endereco das funcoes R para executar o procedimento do ECG Relatorio - Ubuntu
source(file="/home/du/Dropbox/PosDoc/Projetos paralelos/EMG e Torque 2/FuncoesForca.R")
############################################################
#
#     Informacoes descritivas sobre os dados
#       Usuario deve suprir essas informacoes
#
###########################################################

# Diretorios de entrada dos dados brutos
#  e de armazenamento (saida) dos resultados
#  Ambos devem terminar com o simbolo / 

#Para Ubuntu
caminhoEntradaForce<-"/home/du/Dropbox/PosDoc/Projetos paralelos/EMG e Torque 2/Torque/"
caminhoEntradaEMG<-"/home/du/Dropbox/PosDoc/Projetos paralelos/EMG e Torque 2/EMG/"
caminhoSaida<-"/home/du/Dropbox/PosDoc/Projetos paralelos/EMG e Torque 2/Saida/"

#Para Windows
caminhoEntradaForce<-"/Users/Du/Dropbox/PosDoc/Projetos paralelos/EMG e Torque 2/Torque/"
caminhoEntradaEMG<-"/Users/Du/Dropbox/PosDoc/Projetos paralelos/EMG e Torque 2/EMG/"
caminhoSaida<-"/Users/Du/Dropbox/PosDoc/Projetos paralelos/EMG e Torque 2/Saida/"

############################## MUDAR AQUI ##################################################
# Nome da base de dados matlab que contem o miograma   
nomeArqEntrada<- "s186t1m1.xlsx"
Nome<- file_path_sans_ext(nomeArqEntrada)


#Colocar o arquivo em um data frame
DFTotal <- read_excel(paste(caminhoEntradaForce, nomeArqEntrada, sep="/"), col_names=FALSE)
#Arrumando os arquivos e transformando para Newtons.m
DFTorqueN<-separaVetores(DFTotal)
#DFTorqueNLimpo <- DFTorqueN[,"Torque"]


## Colocando o EMG em um data frame
nomeArqEntradaEMG<-paste(Nome,".txt",sep="" )
DFEMG <-read.table(paste(caminhoEntradaEMG, nomeArqEntradaEMG, sep="/"), skip=26,  sep=",")


#Filtrando e passando o RMS
DFEMGFilt<-FiltraRMSJunta(DFEMG)


#Alinhando os vetores
PosicaoTrigger<-min(which(DFTorqueN[, "Trigger"] > max(DFTorqueN[, "Trigger"])*0.8))# Limiar 80%

DFTorqueTemp<-DFTorqueN[PosicaoTrigger:length(DFTorqueN[,"Torque"]),]

CorteET<-min(c(length(DFTorqueTemp[,"Torque"]),length(DFEMGFilt$EMG1)))

DFEMGTorque<-cbind(DFEMGFilt[1:CorteET,], Torque=DFTorqueTemp[1:CorteET,"Torque"])



#Obtendo as informacoes necessarias
# Linha de base
LBase<- mean(DFEMGTorque[1:1000,"Torque"])
#Inicio da contracao
PosicaoIC<-min(which(DFEMGTorque[, "Torque"] > LBase+7.5))
PosicaoFC<-max(which(DFEMGTorque[, "Torque"] > LBase+7.5))

ValorIC<-DFEMGTorque[PosicaoIC, "Torque"]
ValorFC<-DFEMGTorque[PosicaoFC, "Torque"]

#Pico da contracao
PosicaoPicoTorque<-min(which(DFEMGTorque[, "Torque"] == max(DFEMGTorque[, "Torque"])))
PicoContracao<-max(DFEMGTorque[, "Torque"])
PicoContracaoPosIC<- PicoContracao-ValorIC

#Calculo da TDF
#Para Janela 20
Janela<-20
Inclinacoes20<-CalculaTDF(DFEMGTorque,PosicaoIC,Janela,200)

TDF20<-max(Inclinacoes20)
MediaTempoTDF20<-PosicaoIC+Janela*min(which(Inclinacoes20==max(Inclinacoes20)))-(Janela/2)

#Para Janela 1
Janela<-1
Inclinacoes1<-CalculaTDF(DFEMGTorque,PosicaoIC,Janela,200)

TDF1<-max(Inclinacoes1)
MediaTempoTDF1<-PosicaoIC+Janela*min(which(Inclinacoes1==max(Inclinacoes1)))-(Janela)


#Para o calculo do Retardo Eletromecanico
LBaseEMG1<-mean(DFEMGTorque[1000:1500,"EMG1"])
PosicaoICEMG1<-min(which(DFEMGTorque[1000:CorteET, "EMG1"] > LBaseEMG1+0.015))+1000
RetardoEMG1<-PosicaoIC-PosicaoICEMG1

LBaseEMG2<-mean(DFEMGTorque[1000:1500,"EMG2"])
PosicaoICEMG2<-min(which(DFEMGTorque[1000:CorteET, "EMG2"] > LBaseEMG2+0.015))+1000
RetardoEMG2<-PosicaoIC-PosicaoICEMG2

LBaseEMG3<-mean(DFEMGTorque[1000:1500,"EMG3"])
PosicaoICEMG3<-min(which(DFEMGTorque[1000:CorteET, "EMG3"] > LBaseEMG3+0.015))+1000
RetardoEMG3<-PosicaoIC-PosicaoICEMG3


#Integral numerica do EMG 

IntEMG1<-trapz(DFEMGTorque[PosicaoICEMG1:PosicaoFC,"Tempo"],DFEMGTorque[PosicaoICEMG1:PosicaoFC,"EMG1"] ) 
IntEMG2<-trapz(DFEMGTorque[PosicaoICEMG2:PosicaoFC,"Tempo"],DFEMGTorque[PosicaoICEMG2:PosicaoFC,"EMG2"] ) 
IntEMG3<-trapz(DFEMGTorque[PosicaoICEMG3:PosicaoFC,"Tempo"],DFEMGTorque[PosicaoICEMG3:PosicaoFC,"EMG3"] ) 

#Integral e inclinacoes para as diferentes janelas
K<-30
Inclinacoes0_30<-CalculaTDF(DFEMGTorque,PosicaoIC,K,K)
Int_0_30EMG1<-trapz(DFEMGTorque[PosicaoICEMG1:(PosicaoIC+K),"Tempo"], DFEMGTorque[PosicaoICEMG1:(PosicaoIC+K),"EMG1"])
Int_0_30EMG2<-trapz(DFEMGTorque[PosicaoICEMG2:(PosicaoIC+K),"Tempo"], DFEMGTorque[PosicaoICEMG2:(PosicaoIC+K),"EMG2"])
Int_0_30EMG3<-trapz(DFEMGTorque[PosicaoICEMG3:(PosicaoIC+K),"Tempo"], DFEMGTorque[PosicaoICEMG3:(PosicaoIC+K),"EMG3"])

#MAVs
MAV_0_30EMG1<-Int_0_30EMG1/(DFEMGTorque[(PosicaoIC+K),"Tempo"]-DFEMGTorque[PosicaoICEMG1,"Tempo"])
MAV_0_30EMG2<-Int_0_30EMG2/(DFEMGTorque[(PosicaoIC+K),"Tempo"]-DFEMGTorque[PosicaoICEMG2,"Tempo"])
MAV_0_30EMG3<-Int_0_30EMG3/(DFEMGTorque[(PosicaoIC+K),"Tempo"]-DFEMGTorque[PosicaoICEMG3,"Tempo"])

K<-50
Inclinacoes0_50<-CalculaTDF(DFEMGTorque,PosicaoIC,K,K)
Int_0_50EMG1<-trapz(DFEMGTorque[PosicaoICEMG1:(PosicaoIC+K),"Tempo"], DFEMGTorque[PosicaoICEMG1:(PosicaoIC+K),"EMG1"])
Int_0_50EMG2<-trapz(DFEMGTorque[PosicaoICEMG2:(PosicaoIC+K),"Tempo"], DFEMGTorque[PosicaoICEMG2:(PosicaoIC+K),"EMG2"])
Int_0_50EMG3<-trapz(DFEMGTorque[PosicaoICEMG3:(PosicaoIC+K),"Tempo"], DFEMGTorque[PosicaoICEMG3:(PosicaoIC+K),"EMG3"])

#MAVs
MAV_0_50EMG1<-Int_0_50EMG1/(DFEMGTorque[(PosicaoIC+K),"Tempo"]-DFEMGTorque[PosicaoICEMG1,"Tempo"])
MAV_0_50EMG2<-Int_0_50EMG2/(DFEMGTorque[(PosicaoIC+K),"Tempo"]-DFEMGTorque[PosicaoICEMG2,"Tempo"])
MAV_0_50EMG3<-Int_0_50EMG3/(DFEMGTorque[(PosicaoIC+K),"Tempo"]-DFEMGTorque[PosicaoICEMG3,"Tempo"])


K<-100
Inclinacoes0_100<-CalculaTDF(DFEMGTorque,PosicaoIC,K,K)
Int_0_100EMG1<-trapz(DFEMGTorque[PosicaoICEMG1:(PosicaoIC+K),"Tempo"], DFEMGTorque[PosicaoICEMG1:(PosicaoIC+K),"EMG1"])
Int_0_100EMG2<-trapz(DFEMGTorque[PosicaoICEMG2:(PosicaoIC+K),"Tempo"], DFEMGTorque[PosicaoICEMG2:(PosicaoIC+K),"EMG2"])
Int_0_100EMG3<-trapz(DFEMGTorque[PosicaoICEMG3:(PosicaoIC+K),"Tempo"], DFEMGTorque[PosicaoICEMG3:(PosicaoIC+K),"EMG3"])

#MAVs
MAV_0_100EMG1<-Int_0_100EMG1/(DFEMGTorque[(PosicaoIC+K),"Tempo"]-DFEMGTorque[PosicaoICEMG1,"Tempo"])
MAV_0_100EMG2<-Int_0_100EMG2/(DFEMGTorque[(PosicaoIC+K),"Tempo"]-DFEMGTorque[PosicaoICEMG2,"Tempo"])
MAV_0_100EMG3<-Int_0_100EMG3/(DFEMGTorque[(PosicaoIC+K),"Tempo"]-DFEMGTorque[PosicaoICEMG3,"Tempo"])


K<-200
Inclinacoes0_200<-CalculaTDF(DFEMGTorque,PosicaoIC,K,K)
Int_0_200EMG1<-trapz(DFEMGTorque[PosicaoICEMG1:(PosicaoIC+K),"Tempo"], DFEMGTorque[PosicaoICEMG1:(PosicaoIC+K),"EMG1"])
Int_0_200EMG2<-trapz(DFEMGTorque[PosicaoICEMG3:(PosicaoIC+K),"Tempo"], DFEMGTorque[PosicaoICEMG3:(PosicaoIC+K),"EMG3"])
Int_0_200EMG3<-trapz(DFEMGTorque[PosicaoICEMG3:(PosicaoIC+K),"Tempo"], DFEMGTorque[PosicaoICEMG3:(PosicaoIC+K),"EMG3"])

#MAVs
MAV_0_200EMG1<-Int_0_200EMG1/(DFEMGTorque[(PosicaoIC+K),"Tempo"]-DFEMGTorque[PosicaoICEMG1,"Tempo"])
MAV_0_200EMG2<-Int_0_200EMG2/(DFEMGTorque[(PosicaoIC+K),"Tempo"]-DFEMGTorque[PosicaoICEMG2,"Tempo"])
MAV_0_200EMG3<-Int_0_200EMG3/(DFEMGTorque[(PosicaoIC+K),"Tempo"]-DFEMGTorque[PosicaoICEMG3,"Tempo"])


#EMG no Pico
EMGPico1<-trapz(DFEMGTorque[(PosicaoPicoTorque-500):(PosicaoPicoTorque+500),"Tempo"], DFEMGTorque[(PosicaoPicoTorque-500):(PosicaoPicoTorque+500),"EMG1"])
EMGPico2<-trapz(DFEMGTorque[(PosicaoPicoTorque-500):(PosicaoPicoTorque+500),"Tempo"], DFEMGTorque[(PosicaoPicoTorque-500):(PosicaoPicoTorque+500),"EMG2"])
EMGPico3<-trapz(DFEMGTorque[(PosicaoPicoTorque-500):(PosicaoPicoTorque+500),"Tempo"], DFEMGTorque[(PosicaoPicoTorque-500):(PosicaoPicoTorque+500),"EMG3"])


#Construindo um DF para criar a figura
Posicoes<-c(PosicaoIC, PosicaoPicoTorque,MediaTempoTDF20,MediaTempoTDF1)/1000
Valores<-c(ValorIC, PicoContracao,DFEMGTorque[MediaTempoTDF20, "Torque"],DFEMGTorque[MediaTempoTDF1, "Torque"] )
DFAnalises<-data.frame(Posicoes, Valores, row.names=c("Inicio Contracao", "Pico Torque", "TDF20", "TDF1"))

Posicoes2<-c(PosicaoIC,PosicaoFC, PosicaoPicoTorque,MediaTempoTDF20,MediaTempoTDF1)/1000
Valores2<-c(ValorIC, ValorFC,PicoContracao,DFEMGTorque[MediaTempoTDF20, "Torque"],DFEMGTorque[MediaTempoTDF1, "Torque"] )
DFAnalises2<-data.frame(Posicoes2, Valores2, row.names=c("Inicio Contracao","Fim da Contracao", "Pico Torque", "TDF20", "TDF1"))

#Para colocar as linhas verticais

segment_data = data.frame(
  x = c(DFEMGTorque[PosicaoICEMG1, "Tempo"],DFEMGTorque[PosicaoIC+30, "Tempo"], DFEMGTorque[PosicaoIC+50, "Tempo"],DFEMGTorque[PosicaoIC+100, "Tempo"],DFEMGTorque[PosicaoIC+200, "Tempo"]),
  xend = c(DFEMGTorque[PosicaoICEMG1, "Tempo"],DFEMGTorque[PosicaoIC+30, "Tempo"], DFEMGTorque[PosicaoIC+50, "Tempo"],DFEMGTorque[PosicaoIC+100, "Tempo"],DFEMGTorque[PosicaoIC+200, "Tempo"]),
  y = c(DFEMGTorque[LBase, "Torque"],DFEMGTorque[LBase, "Torque"],DFEMGTorque[LBase, "Torque"],DFEMGTorque[LBase, "Torque"],DFEMGTorque[LBase, "Torque"]),
  yend = c(DFEMGTorque[PosicaoIC, "Torque"],DFEMGTorque[PosicaoIC+30, "Torque"],DFEMGTorque[PosicaoIC+50, "Torque"],DFEMGTorque[PosicaoIC+100, "Torque"],DFEMGTorque[PosicaoIC+200, "Torque"])
)



# Figura Torque
Fig<- ggplot(data=DFEMGTorque, aes(x=Tempo, y=Torque)) + geom_line( size = 3) + ylab("Torque (N.m)")+ xlab("Tempo (segundos)") 
Fig<- Fig + theme_bw()  +  theme(legend.position="none")
Fig<- Fig + geom_point(data=DFAnalises2, aes(x = Posicoes2, y = Valores2, color=rownames(DFAnalises2)), size=30)
Fig<- Fig + theme(axis.text=element_text(size=30*2.5), axis.title=element_text(size=30*2.5))  + theme( axis.line = element_line(size = 3, linetype = "solid"))
#Fig+ theme(legend.text=element_text(size=50))


FigZoom<- ggplot(data=DFEMGTorque[1500:PosicaoPicoTorque+600,], aes(x=Tempo, y=Torque)) + geom_line( size = 3) + ylab("Torque (N.m)")+ xlab("Tempo (segundos)") 
FigZoom<- FigZoom + theme_bw()  +  theme(legend.title=element_blank())
FigZoom<- FigZoom + geom_point(data=DFAnalises, aes(x = Posicoes, y = Valores, color=rownames(DFAnalises)), size=30)
FigZoom<- FigZoom + geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend), size=3)
FigZoom<- FigZoom + geom_segment(aes(x = DFEMGTorque[PosicaoPicoTorque-500, "Tempo"], y = PicoContracao+5, xend = DFEMGTorque[PosicaoPicoTorque+500, "Tempo"], yend = PicoContracao+5), size=5, color="red")
FigZoom<- FigZoom + theme(axis.text=element_text(size=30*2.5), axis.title=element_text(size=30*2.5),legend.text=element_text(size=50))  + theme( axis.line = element_line(size = 3, linetype = "solid"))
#FigZoom


# Figura EMG1
FigEMG1<- ggplot(data=DFEMGTorque, aes(x=Tempo, y=EMG1)) + geom_line( size = 1.5) + ylab("EMG 1 (mV)")+ xlab("Tempo (segundos)") 
FigEMG1<- FigEMG1 + theme_bw()  +  theme(legend.title=element_blank())
FigEMG1<- FigEMG1 + geom_point(aes(x=DFEMGTorque[PosicaoICEMG1,"Tempo"],y=DFEMGTorque[PosicaoICEMG1,"EMG1"]), color="red", size=20) 
FigEMG1<- FigEMG1 + theme(axis.text=element_text(size=30*2.5), axis.title=element_text(size=30*2.5))  + theme( axis.line = element_line(size = 3, linetype = "solid"))
#FigEMG1

# Figura EMG2
FigEMG2<- ggplot(data=DFEMGTorque, aes(x=Tempo, y=EMG2)) + geom_line( size = 1.5) + ylab("EMG 2 (mV)")+ xlab("Tempo (segundos)") 
FigEMG2<- FigEMG2 + theme_bw()  +  theme(legend.title=element_blank())
FigEMG2<- FigEMG2 + geom_point(aes(x=DFEMGTorque[PosicaoICEMG2,"Tempo"],y=DFEMGTorque[PosicaoICEMG2,"EMG2"]), color="red", size=20) 
FigEMG2<- FigEMG2 + theme(axis.text=element_text(size=30*2.5), axis.title=element_text(size=30*2.5))  + theme( axis.line = element_line(size = 3, linetype = "solid"))
#FigEMG2

# Figura EMG3
FigEMG3<- ggplot(data=DFEMGTorque, aes(x=Tempo, y=EMG3)) + geom_line( size = 1.5) + ylab("EMG 3 (mV)")+ xlab("Tempo (segundos)") 
FigEMG3<- FigEMG3 + theme_bw()  +  theme(legend.title=element_blank())
FigEMG3<- FigEMG3 + geom_point(aes(x=DFEMGTorque[PosicaoICEMG3,"Tempo"],y=DFEMGTorque[PosicaoICEMG3,"EMG3"]), color="red", size=20) 
FigEMG3<- FigEMG3 + theme(axis.text=element_text(size=30*2.5), axis.title=element_text(size=30*2.5))  + theme( axis.line = element_line(size = 3, linetype = "solid"))
#FigEMG3

png(file=paste(caminhoSaida,Nome,".png",sep=""), width=2220*2, height=1850*2.5)
grid.arrange(FigZoom, Fig, FigEMG1, FigEMG2, FigEMG3, ncol=1)
dev.off()



#Colocando os resultados em uma tabela 
tabela<-data.frame(Nome, ValorIC, PicoContracao, PicoContracaoPosIC,TDF20,TDF1,RetardoEMG1,RetardoEMG2,RetardoEMG3,IntEMG1,IntEMG2,IntEMG3, Inclinacoes0_30,Inclinacoes0_50,Inclinacoes0_100,Inclinacoes0_200,Int_0_30EMG1,Int_0_30EMG2,Int_0_30EMG3,Int_0_50EMG1,Int_0_50EMG2,Int_0_50EMG3,Int_0_100EMG1,Int_0_100EMG2,Int_0_100EMG3,Int_0_200EMG1,Int_0_200EMG2,Int_0_200EMG3,Inclinacoes0_30,MAV_0_30EMG1,MAV_0_30EMG2,MAV_0_30EMG3,MAV_0_50EMG1,MAV_0_50EMG2,MAV_0_50EMG3,MAV_0_100EMG1,MAV_0_100EMG2,MAV_0_100EMG3,MAV_0_200EMG1,MAV_0_200EMG2,MAV_0_200EMG3,EMGPico1,EMGPico2,EMGPico3)
#Para o primeiro
#write.table(tabela,paste(caminhoSaida,"Resultados", ".csv", sep=""),sep=",", row.names = FALSE)

#Para os outros
write.table(tabela,paste(caminhoSaida,"Resultados", ".csv", sep=""),append=T, sep=",", row.names = FALSE, col.names = FALSE)


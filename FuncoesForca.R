# NOME: separaVetores
# DESC: Separa um DF vindo do excel em um data frame com os dados de torque (ja transformado em Newtons.m)
# PARAMETROS:
#   Dataframe [DF]: Data frame proveniente do excel
# EXEMPLO DE CHAMADA:
#   resultado <- separaVetores(Dataframe)

separaVetores<- function(DFTotal){

  #Separando em vetores de torque e dados para transformar em Newtons.
  DFTorque<-DFTotal[seq(1,5)]
  DFTransform<-DFTotal[seq(7,13)]
  colnames(DFTransform) <- DFTransform[1,]
  DFTransform <- na.omit(DFTransform[-1, ] )
  DFTransform <- as.data.frame(sapply(DFTransform, as.numeric)) 
  
  #transformando em Newton
  DFTempN<-(DFTorque[ ,"X__2"]-DFTransform["peso da perna",])*(DFTransform["Peso em N.m", 1]/abs(DFTransform["calibração",]-DFTransform["sem nada",]))
  #head(DFTorqueN)
  DFTorqueN<-DFTorque
  DFTorqueN[ ,"X__2"]<-DFTempN
  
  #Adicionando um vetor de tempo e nomes para colunas
  DFTorqueN<-cbind(T=seq(1/1000, nrow(DFTorqueN)/1000, by=1/1000), DFTorqueN)
  colnames(DFTorqueN) <- c("Tempo", "Trigger1","Torque","EMG","TriggerUS","Trigger")
  
  return(DFTorqueN)
}

# NOME: CalculaTDF
# DESC: Calcula a taxa de desenvolvimento de Forca para um data frame de Torque ao longo do tempo
# PARAMETROS:
#   Dataframe [DF]: Data frame de Torque ao longo do tempo
#   PosicaoIC [INT]: Posicao do inicio da contracao
#   Janela [INT]: Janela para calcular a inclinacao
#   Total [INT]: Numero de pontos total onde a Janela sera passada, se Janela==Tempo o resultado será apenas uma inclinacao.
# EXEMPLO DE CHAMADA:
#   resu0ltado <- CalculaTDF(Dataframe,PosicaoIC, Janela, Total )

CalculaTDF<- function(DFTorqueN, PosicaoIC, Janela, Total){
  j<-1
  Inclinacoes<-rep(0,Total/Janela)
  for (i in seq(PosicaoIC,PosicaoIC+(Total-Janela),Janela)){
    Inclinacoes[j]<-(DFTorqueN[i+Janela, "Torque"]-DFTorqueN[i, "Torque"])/(DFTorqueN[i+Janela, "Tempo"]-DFTorqueN[i, "Tempo"])
    j<-j+1
  }
  return(Inclinacoes)
}


# NOME: FiltraRMS
# DESC: Filtra o vetor utilizando um filtro Butterworth com ordem 4, zero-lag e frequencia de corte de 5Hz. Em sequencia passa um rms com janela 50
# PARAMETROS:
#   Vetor [num]: Vetor de EMG ao longo do tempo
#
# EXEMPLO DE CHAMADA:
#   resultado <- FiltraRMS(Vetor)

FiltraRMS<- function(DF_EMG_separado){
  
  DF_EMGL<-subset(DF_EMG_separado, DF_EMG_separado<2)
  DF_EMGL<-DF_EMGL*1000 #Para transformar volts em milivolts
  bf <- butter(4, 5/500, type = "high")
  b <- filter(bf, DF_EMGL)
  bzlrev<-filter(bf, rev(b))
  bzl<-rev(bzlrev)
  
  j<-1
  emgrms<-rep(0,(length(bzl)-49))
  for (i in seq(1,length(bzl)-49,1)){
    emgrms[j]<-sqrt(mean(bzl[i:i+49]^2))
    j<-j+1
  }
  return(emgrms)
}

# NOME: FiltraRMSJunta
# DESC: Filtra o vetor utilizando um filtro Butterworth com ordem 4, zero-lag e frequencia de corte de 5Hz. Em sequencia passa um rms com janela 50 e junta tudo em um Data Frame
# PARAMETROS:
#    Vetor [num]: Vetor de EMG ao longo do tempo
# EXEMPLO DE CHAMADA:
#   resultado <- FiltraRMSJunta(Dataframe)

FiltraRMSJunta<- function(DFEMG_import){
  EMG_1<- FiltraRMS(DFEMG_import[,2])
  EMG_2<- FiltraRMS(DFEMG_import[,3])
  EMG_3<- FiltraRMS(DFEMG_import[,4])
  Tempo<-DFEMG_import[, 1]
  Corte<-min(c(length(Tempo),length(EMG_1),length(EMG_2),length(EMG_3)))
  
  DFEMGFiltrado<- data.frame(Tempo[1:Corte],EMG_1[1:Corte], EMG_2[1:Corte], EMG_3[1:Corte])
  #Adicionando nomes para colunas
  colnames(DFEMGFiltrado) <- c("Tempo", "EMG1","EMG2","EMG3")
  
  return(DFEMGFiltrado)
}


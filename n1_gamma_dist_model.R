## Hacemos una distribución gamma a partir de TODAS las muestras que tengamos a 
# partir de t=8. Consideramos que TODOS LOS TIEMPOS SON EQUIVALENTES

# Neutral model (Goyal paper)

## GUARDO UNA TABLA PARA CADA TIME POINT PARA CADA M_INIC
################################################################################
# Real abundances follow a gamma distribution. So we're going to use that
# to infer simulated abundances.
# - FOR EACH SPECIES for each timepoint, take the mean abundance from the real data
# - and the variance
library("gsubfn") # for unpacking values
source("./my_functions.R")

exp_f="./data/all_transfers_table_glc.txt"
map_f="./data/map_glc_ALL.csv"

exp<-get_abundance_and_tax_from_table(exp_f)[[1]] # my_functions.R
# Filtro: solo simularé aquellas OTUs que no se extinguen en los valores reales   # !!!
exp <- exp[exp[map[map["T"]==12,][,"SA"]] %>% rowSums()!=0,]                   

# need the map for grouping by m_inic
map <- read.csv(map_f,sep=",")
# solo podemos usar X2 y X6
m_inic <- unique(map["ORIG"])[["ORIG"]]

# simulation parameters
n = 1000          # number of simulations                                 
timepoints = 5    # number of time points                                   

# ==============================================================================
#   Per m_inic
# ==============================================================================
simuls <- list()
for (m in m_inic){
  # Vamos a hacer una distribución gamma a partir de todos los timepoints estables
  sel=map[map["ORIG"]==m & map["T"]>7,]
  exp_m <- exp[sel[,"SA"]]

  for (t in 1:timepoints) {
    # gamma simulations for all the time points
    simuls[[m]][[t]] <-
      create_gamma_distr_simuls(data=exp_m, n=n) %>% as.data.frame() %>% my_transpose()
  }
}

# 4. Save data
for (m in m_inic){
  for (t in 1:timepoints) {
    # save (a file per time point)
    write.csv(x = simuls[[m]][[t]] %>% my_transpose,file = paste0("n1_results","/simul_",m,"_t",t,".csv"))
  }
}

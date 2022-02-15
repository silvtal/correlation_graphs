# - For the null model, a gamma distribution is constructed for each OTU based on the experimental data. 
# 
# - Then, random communities are generated from the gamma distributions and community compositions are renormalized by dividing each individual abundance by the communities' total sum. 
# 
# - The simulated communities are then arranged in passage trajectories and the correlations between each pair of OTUs obtained as per the experimental data. 
# 
# - Finally, the P values of the experimentally observed correlations are obtained by comparing with the expected null distribution of 1000 simulated correlations. 

# 13feb2022 X2 and X6 have multiple replicates. If only_one_rep==TRUE, we'll 
#           simply choose one. If FALSE, we will obtain correlation data for
#           each replicate separately.

# ------------------------------------------------------------------------------
# 13feb2022 Para X2 y X6 vamos a usar solo una de las muchas réplicas, para
#           ser coherentes / que sea más fácil de explicar. Así que si 
#           only_one_rep == TRUE, guardará solo la primera que haya

# WARNING : no confundir este script con el n2_correlations.R de mi null model !!
# el que quiero es, en realidad, ESTE!
# ------------------------------------------------------------------------------


# Libraries, parameters, wd
# ==============================================================================
setwd("/home/silvia/AAA/2021-09-22_gamma_model_networks/")

library("gsubfn") # for unpacking values
source("/home/silvia/Apps/functions_for_neutral_modelling.R") # source scripts
source("/home/silvia/Apps/my_functions.R")                    # are in neutral_model
library("plyr")

exp_f="/home/silvia/AAA/OTU_data/all_transfers_table_glc.txt"
map_f="/home/silvia/AAA/data/map_glc_ALL.csv"
simul_input_file = "gamma_model_glc"
simul_folder = "./gamma_simul_correlations"
real_folder  = "./gamma_real_correlations"

m_inic = paste0("X",c(1:12))

n = 1000            # number of simulations                                        # !!! important to get right.
timepoints_S = 1:5  # time points for simulations

reps = 1:8          # number of replicates (no incluir la 0 !)
timepoints_R = 8:12 # time points for real data
                    # Nos interesan las comunidades estables
only_one_rep = TRUE




# Cargar datos reales
# ==============================================================================
exp<-get_abundance_and_tax_from_table(exp_f)[[1]] # my_functions.R
map <- read.csv(map_f,sep=",")

## Quité este filtro en pos del filtro de "mínimo 10"
# # Filtro 1: solo miro las correlaciones entre aquellas OTUs que no se extinguen
# # (que no sean 0 en todas las transfers número 12)
# exp <- exp[exp[map[map["T"]==12,][,"SA"]] %>% rowSums()!=0,] %>% my_transpose()                

# Creamos all_reps, que es equivalente a model_simul. Es una lista de m_inic,     # quedan fuera 4, 5 y 9 !!!!! les falta la transfer 11 de la replica 4 (que es la que supuestamente tiene siempre los tr de 8 a 12)
# siendo cada m_inic una lista de 5 time points, siendo cada time point un DF 
# con tantas filas como réplicas y tantas columnas como OTUs
all_reps <- list()  
all_kept_otus <- list()
for (m in m_inic){
  at_least_one_good_rep = FALSE
  for (rep in reps) {
    selected_sa = map[map["ORIG"] == m &                    # in that m_inic
                        map["REP"]== rep &                  # in the given rep
                        map["T"] >=min(timepoints_R),][,"SA"] # select only stable   
    
    ## Importante seleccionar réplicas con 5 time points para poder sacar una 
    ## correlación significativa (es casi siempre la replica 4, salvo para las 
    ## replicas de las muestras iniciales X2 y X6, que todas tienen los 5)
    if (length(selected_sa)==5) {
      # Filtro: solo miro EN CADA REP/SIMUL las correlaciones entre aquellas OTUs
      # que en al menos 3 de los timepoints a estudiar superan 10 reads   ## FIXED nov 8-12
      # exp <- exp[exp[map[map["T"]==12,][,"SA"]] %>% rowSums()!=0,] %>% my_transpose()                
      # exp <- exp[exp[map[map["T"]>=min(timepoints_R),][,"SA"]] != 0,] %>% my_transpose()
      selected_exp <- exp[selected_sa]# / colSums(exp[selected_sa]) # abs abu -> RA
      keep <- apply(selected_exp, 1, function(r) {
        th <- 10 # threshold
        return(sum(r > th)>3)
      })
      selected_exp <- selected_exp[names(keep[keep]),] # abs abu
      if (length(selected_sa)>0){
        all_kept_otus[[m]] <- c(all_kept_otus[[m]], rownames(selected_exp)) # for later
      }
      all_reps[[m]][[as.character(rep)]] <- selected_exp %>% my_transpose # guardo los datos de esa rep
      found_at_least_one_good_rep = TRUE
    }
  }
  if (!found_at_least_one_good_rep) {message(paste("Didn't find any replicates with 5 registered transfers/timepoints for sample",m))}
}


# Medir correlaciones reales y guardarlas
# ==============================================================================
corrs = list()
message("Obteniendo valores de correlación y p-valores (REALES)")
for (m in names(all_reps)){ # no son todas las m_inic !
  message(paste0("Analizando réplicas de la muestra inicial: ",m))
  for (rep in names(all_reps[[m]])) {
    rep = as.character(rep) #~
    # Las réplicas aquí correrán el mismo destino que las 
    # simulaciones de los datos simulados. El DF final contiene todos los time 
    # points de esa réplica y a partir de ellos se harán las correlaciones.
    message(paste0("Réplica ",rep))
    # rep_tp == DF con todas las timepoints de la réplica "rep"
    rep_tp <- all_reps[[m]][[as.character(rep)]]
    
    if (nrow(rep_tp) > 1) { # solo puedo sacar pares si hay al menos dos OTUs
      message("Computing Pearson correlations...")
      pairs = t(combn(colnames(rep_tp),2)) %>% as.data.frame()
      # corrs == list. Dentro de corrs$m$PCG$rep tendremos dos DF, uno con 
      #                p-vals y otro con las corrs en sí. Tienen tantas cols y
      #                rows como OTUs
      
      # DF vacíos:
      corrs[[m]][[rep]]<-list()
      corrs[[m]][[rep]]$rep_pvals <- data.frame(row.names = colnames(rep_tp))
      corrs[[m]][[rep]]$rep_corrs <- data.frame(row.names = colnames(rep_tp))
      
      # Hago Pearson para cada pareja y meto los resultados en las casillas 
      # correspondientes de los DFs
      for (row in 1:nrow(pairs)) {
        p=pairs[row,]
        s=cor.test(rep_tp[,p[[1]]], rep_tp[,p[[2]]], method=c("pearson"))          # quiza no lo mas eficiente pero no importa...
        corrs[[m]][[rep]]$rep_pvals[p[[1]],p[[2]]] <- s$p.value
        corrs[[m]][[rep]]$rep_pvals[p[[2]],p[[1]]] <- s$p.value
        corrs[[m]][[rep]]$rep_corrs[p[[1]],p[[2]]] <- s$estimate[["cor"]]
        corrs[[m]][[rep]]$rep_corrs[p[[2]],p[[1]]] <- s$estimate[["cor"]]
      }
      # Order the output...
      A <- corrs[[m]][[rep]]$rep_corrs
      corrs[[m]][[rep]]$rep_corrs <- A[attr(A, 'row.names'), attr(A, 'row.names')]
      A <- corrs[[m]][[rep]]$rep_pvals
      corrs[[m]][[rep]]$rep_pvals <- A[order(as.character(attr(A, 'row.names'))), order(colnames(A))]

    } else {
      message("Less than 2 OTUs pass the filter. Skipping to the next replicate...")
    }
  }
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # FILTRO;nov10: como aplico un filtro, tenemos diferentes otus seleccionadas   # TODO está feo programado... aunque va perfecto
  # en cada matriz/replica. Tengo que rellenar con ceros para que todas tengan
  # las mismas otus en las row/colnames
  for (rep in names(corrs[[m]])) { # solo las que tienen más de 1 OTU !
    kept_otus <- all_kept_otus[[m]]%>%unique # todas las otus seleccionadas en
                                             # todas las reps de esta m_inic
    temp <- setNames(data.frame(matrix(ncol=length(kept_otus),nrow=0)),
                     kept_otus) # creo df vacío
    
    temp[kept_otus,] <- corrs[[m]][[rep]]$rep_corrs[kept_otus,] # meto datos
    temp[is.na(temp)] <- 0 # meto los ceros en los huecos que quedan
    corrs[[m]][[rep]]$rep_corrs <- temp
    # same for pvals
    temp[kept_otus,] <- corrs[[m]][[rep]]$rep_pvals[kept_otus,] # meto datos
    temp[is.na(temp)] <- 1 # meto los UNOS  en los huecos que quedan
    corrs[[m]][[rep]]$rep_pvals <- temp
  } # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  

  # Cuando tenemos los dos DF llenos para cada rep, tengo que crear dos 
  # nuevos DF comunes a todas las reps, que contengan los valores medios de
  # correlacion y de p-values
  # (Si no se ha podido llenar porque no hay correlaciones válidas en ninguna de
  #  las trayectorias de esa rep, debido a que no hay suficientes OTUs pasando
  #  el filtro, no hay nada que guardar !)
  if (length(corrs[[m]])!=0){
    # Guardo una matriz de corrs y otra de pvalor para cada m_inic
    all_corrs = lapply(corrs[[m]], function(rep){rep$rep_corrs})    
    all_pvals = lapply(corrs[[m]], function(rep){rep$rep_pvals})
    
    if (length(all_corrs)>1) { # solo puedo hacer medias si hay más de dos matrices !
      # 13feb2022 Para X2 y X6 vamos a usar solo una de las muchas réplicas, para
      #           ser coherentes / que sea más fácil de explicar. Así que si 
      #           only_one_rep == TRUE, guardará solo la primera que haya
      if (!only_one_rep) {
        aaply(laply(.data = all_corrs, as.matrix), 
              .margins = c(2,3), 
              .fun = mean) %>%
          write.csv(paste0(real_folder,"/",m,"_all_corrs.csv"))
        aaply(laply(.data = all_pvals, as.matrix), 
              .margins = c(2, 3), 
              .fun = mean) %>%
          write.csv(paste0(real_folder,"/",m,"_all_pvals.csv"))
      } else {
        all_corrs[[names(all_corrs)[1]]] %>%
          write.csv(paste0(real_folder,"/",m,"_all_corrs.csv"))
        all_pvals[[names(all_pvals)[1]]] %>%
          write.csv(paste0(real_folder,"/",m,"_all_pvals.csv"))
      }
    } else {
      all_corrs[[names(all_corrs)]] %>%
        write.csv(paste0(real_folder,"/",m,"_all_corrs.csv"))
      all_pvals[[names(all_pvals)]] %>%
        write.csv(paste0(real_folder,"/",m,"_all_pvals.csv"))
    }
    # Limpio memoria
    corrs[[m]] = NULL
  } else {
    message(paste0("Could not compute correlations for ", m))
  }
}

# Load all gamma simulated data (there's a dataframe for each time step)
# ==============================================================================
gamma_simul = list()
for (m in m_inic){
  for (t in timepoints_S) {
    t=as.character(t)
    gamma_simul[[m]][[t]]<-read.csv(paste0(simul_input_file,"/simul_",m,"_t",t,".csv"),
                                    row.names = 1, check.names = F)
  }
}

# Then, for each inter-species pair, we then measured the Pearson correlation 
# coefficient between their abundance trajectories in each time point
# - Todos los time points de todas las réplicas (reales) y de todas las simuls,
#   pero separando a cada réplica/iteración.
# - Luego hago la media de las correlaciones obtenidas.
# - Reales por un lado y simulaciones por otro lado. Solo tiempos finales.
# ==============================================================================
corrs = list()
all_kept_otus_s <- list()
message("Obteniendo valores de correlación y p-valores (SIMULACIONES)")
for (m in names(all_reps)){ # no son todas las m_inic                            # !!! Nos ahorramos calcular correlaciones para m_inic cuyos datos reales no tenemos. renta, creo, pero si me da una millonésima parte de problema lo quito
  message(paste0("Analizando simulaciones a partir de réplicas de la muestra inicial: ",m))
  for (simul in 1:n) {# Ahora para cada simulacion voy leyendo su fila 
                      # correspondiente para todos los time points. El DF
                      # final contiene todas las time points de esa simula-
                      # ción, y a partir de ellos se harán las correlaciones.
                      # Finalmente se hace la media de las correlaciones de
                      # todas las simulaciones.
    simul = as.character(simul)
    
    message(paste0("Simulación ",simul))
    # g == DF con todas las timepoints de la simulación "simul"
    g <- data.frame(row.names=timepoints_S) # la vaciamos en cada simul
    for (t in timepoints_S) {
      t = as.character(t)
      g <- rbind(g, gamma_simul[[m]][[t]][simul,])
    }
    
    # Filtro: PARA CADA RÉPLICA/SIMUL, solo miro las correlaciones entre aquellas
    # OTUs que en al menos 3 de los 5 timepoints a estudiar superan 10 reads     ## FIXED nov 8-12 + added to simuls jan 20 2022
    keep <- apply(g, 2, function(r) {
      th <- 10 # threshold
      return(sum(r > th)>3)
    })
    g <- g[names(keep[keep])] # abs abu
    all_kept_otus_s[[m]] <- c(all_kept_otus_s[[m]], colnames(g)) # for later
    
    # solo si más de 1 OTU han pasado el filtro podemos continuar y hacer pares
    # para esa simulación
    if (ncol(g) > 1) {
      pairs = t(combn(colnames(g),2)) %>% as.data.frame()
      # corrs == list. Dentro de corrs$m$PCG$simul tendremos dos DF, uno con 
      #                p-vals y otro con las corrs en sí. Tienen tantas cols y
      #                rows como OTUs
    
      # DF vacíos:
      corrs[[m]][[simul]]<-list()
      corrs[[m]][[simul]]$simul_pvals <- data.frame(row.names = colnames(g))
      corrs[[m]][[simul]]$simul_corrs <- data.frame(row.names = colnames(g))
      
      # Hago Pearson para cada pareja y meto los resultados en las casillas 
      # correspondientes de los DFs
      message("Computing Pearson correlations...")
      for (row in rownames(pairs)) {
        p=pairs[row,]
        s=cor.test(g[,p[[1]]], g[,p[[2]]], method=c("pearson"))
        corrs[[m]][[simul]]$simul_pvals[p[[1]],p[[2]]] <- s$p.value
        corrs[[m]][[simul]]$simul_pvals[p[[2]],p[[1]]] <- s$p.value
        corrs[[m]][[simul]]$simul_corrs[p[[1]],p[[2]]] <- s$estimate[["cor"]]
        corrs[[m]][[simul]]$simul_corrs[p[[2]],p[[1]]] <- s$estimate[["cor"]]
      }
      # Order the output...
      A <- corrs[[m]][[simul]]$simul_corrs
      corrs[[m]][[simul]]$simul_corrs <- A[attr(A, 'row.names'), attr(A, 'row.names')]
      A <- corrs[[m]][[simul]]$simul_pvals
      corrs[[m]][[simul]]$simul_pvals <- A[order(as.character(attr(A, 'row.names'))), order(colnames(A))]
    } else {
      message("Less than 2 OTUs pass the filter. Skipping to the next simulation...")
    }
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # FILTRO;nov10: como aplico un filtro, tenemos diferentes otus seleccionadas   # TODO está feo programado... aunque va perfecto
  # en cada matriz/simul Tengo que rellenar con ceros para que todas tengan
  # las mismas otus en las row/colnames
  for (simul in names(corrs[[m]])) { # solo las que tienen más de 1 OTU !
    kept_otus_s <- all_kept_otus_s[[m]]%>%unique # todas las otus seleccionadas en
    # todas las simuls de esta m_inic
    temp <- setNames(data.frame(matrix(ncol=length(kept_otus_s),nrow=0)),
                     kept_otus_s) # creo df vacío
    
    temp[kept_otus_s,] <- corrs[[m]][[simul]]$simul_corrs[kept_otus_s,] # meto datos
    temp[is.na(temp)] <- 0 # meto los ceros en los huecos que quedan
    corrs[[m]][[simul]]$simul_corrs <- temp
    # same for pvals
    temp[kept_otus_s,] <- corrs[[m]][[simul]]$simul_pvals[kept_otus_s,] # meto datos
    temp[is.na(temp)] <- 1 # meto los UNOS  en los huecos que quedan
    corrs[[m]][[simul]]$simul_pvals <- temp
  } # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    # Cuando tenemos los dos DF llenos para cada simul, tengo que crear dos        # FIXED 3oct2021
    # nuevos DF comunes a todas las simuls, que contienen los 1000 valores de
    # correlacion y de p-values simulados para cada pareja posible de OTUs.
    # (En el siguiente script, n3, sacaré una distribución para cada pareja)
    
    # Guardo una matriz de 1000 corrs y otra de 1000 pvalor por OTU pair para 
    # cada m_inic
    # Primero una lista de matrices
    all_corrs = lapply(corrs[[m]], function(simul){simul$simul_corrs})
    all_pvals = lapply(corrs[[m]], function(simul){simul$simul_pvals})
    
    # Y luego junto y simplifico
    flattenCorrMatrixCor(all_corrs) %>% 
      write.csv(paste0(simul_folder,"/",m,"_all_corrs.csv"))
    
    flat_pvals <- flattenCorrMatrixCor(all_pvals) %>%
      write.csv(paste0(simul_folder,"/",m,"_all_pvals.csv"))
    
    # Limpio memoria
    corrs[[m]] = NULL
}
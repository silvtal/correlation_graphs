# The "real correlations" here are the means of the correlations obtained
# from a selection of replicates+timesteps. We selected only those replicates
# which had the last five transfers (=stable community). For each replicate,
# we made a correlation matrix from data from those 5 timesteps. Our "real
# correlations" are obtained as the mean from the correlation matrices of all
# available replicates.
## ^ unless only_one_rep==TRUE during n2_correlations.R

# Aditionally, for the graphs we'll select only those correlations that are 
# significant according to a p-value obtained through a permutation test.
# ------------------------------------------------------------------------------

# Las "correlaciones reales" que analizaremos son las correlaciones medias 
# hechas a partir de una selección de réplicas y tiempos. Solo se 
# seleccionaron aquellas réplicas con los 5 últimos transfers. Para cada 
# réplica se hizo una matriz de correlación a partir de esos 5 transfers 
# (time points) y luego la media de estas matrices de correlación es la
# que cargamos.

# A su vez, vamos a seleccionar solo aquellas correlaciones "reales" que sean 
# significativas de acuerdo a un p-valor de pertenencia a la distribución de 
# correlaciones simuladas a partir de la gamma

# Por último graficamos y guardamos grafos de correlación

# ------------------------------------------------------------------------------

# Libraries, parameters, wd
# ==============================================================================
setwd("/home/silvia/AAA/2021-09-22_gamma_model_networks/")

library("gsubfn") # for unpacking values
source("/home/silvia/Apps/functions_for_neutral_modelling.R")
source("/home/silvia/Apps/my_functions.R")
if("tidyverse" %in% (.packages())){
  detach("package:tidyverse", unload=TRUE)
}
library("plyr")
library("tidyverse")
# visualize corrs
library("corrplot")
library("igraph")
library("gridGraphics")
library("grid")
library("ggplotify")
library("ggpubr")

simul_folder = "./gamma_simul_correlations"
real_folder  = "./gamma_real_correlations"
outputdir = "./gamma_model_GRAPHS"

taxaf = "../OTU_data/all_transfers_table_glc.txt"
taxa<-get_abundance_and_tax_from_table(taxaf)[[2]] 

m_inic = paste0("X",1:12) # irá probando· Si no hay para esa muestra inicial, 
                          # se dará cuenta. Lo mismo para las réplicas           # fixed 11 oct 
mypval = 0.025
# For each m_inic
# ==============================================================================
for (m in m_inic) {
  tryCatch(
    expr={
      ## Load all data
      ## =============
      real_corrs <- read.csv(paste0(real_folder,"/",m,"_all_corrs.csv"),check.names=FALSE,row.names = 1) %>% 
        replace(is.na(.),0) %>% flattenCorrMatrix
      
      simul_corrs <- read.table(paste0(simul_folder,"/",m,"_all_corrs.csv"),
                                check.names=FALSE,row.names = 1,sep = ",",
                                header = T) %>% 
        replace(is.na(.),0) 
      
      # simul_pvals <- read.table(paste0(simul_folder,"/",m,"_all_pvals.csv"),
      #                           check.names=FALSE,row.names = 1,sep = ",",
      #                           header = T) %>%
      #   replace(is.na(.),0)
      
      ## Check p value of real correlation according to simul distribution
      ## =================================================================
      ## !!! importante: esto se puede hacer así solo si las primeras dos columnas
      ##     de la matriz real son iguales que las de la matriz simul
      n = ncol(simul_corrs)
          ## Permutation test: you take your "original" data and sample as many 
          ## values as you have in your "experimental sample". We have one (we test
          ## each real correlation separately). We do this 1000 times.
          # sampled <- c()
          # for (i in 1:1000) {
          #   sampled[i] <- sample(simul_corrs[c,3:n]%>%as.numeric, size = 1)
          # }
          ##     --> we've actually already done this : our simul_corrs are a 1
          ##         value sampling, repeated 800 or 1000 times !!!
          ## So, for each real value: 
            # How often is the real value bigger than the sampled one?
            # How often is it lower?
            # ^if any of these two values is less than 0.025, we assume it's not in
            #  the null distribution.      
      
      for (c in 1:nrow(real_corrs)){
        perms   <- simul_corrs[c,3:n]%>%as.numeric
        ## to real value to check
        x    <- real_corrs[c,3]
        ## obtain and save p-val
        pval <- min(sum(x<=perms), sum(x>=perms))/length(perms)
        # print(paste("sum(x<=perms)", sum(x<=perms))) # DEBUG
        # print(paste("sum(x<=perms)", sum(x>=perms)))
        real_corrs[c,"pvalue"] <- pval
      }
  
      real_corrs["pvalue"][is.na(real_corrs["pvalue"])]<-1
      
      # Cambio colnames
      colnames(real_corrs) <- c("otu1","otu2","weight","pvalue")

      # Quito las parejas con correlación 0
      real_corrs<-real_corrs[real_corrs$weight!=0,]
      
      # (plot2) plot the "corrplot" 
      ## ===========================================================================
      df=unflatten(real_corrs) #functions_for_neutral_modelling.R
      c_scale = c('#7F0000', 'red', '#FF7F00', 'yellow', 'white',
                  'cyan', '#007FFF', 'blue', '#00007F')%>% rev
      col2 = colorRampPalette(c_scale)
      # library(RColorBrewer); coln = brewer.pal(n = 10, name = 'PRGn')
      corrplot(corr = df%>%as.matrix(), type = 'lower', order = 'hclust', tl.col = 'black', 
               cl.ratio = 0.2, tl.srt = 45, col = col2(10), bg = "lightgrey",
               sig.level = mypval, p.mat = unflatten(real_corrs,values="pvalue")%>%as.matrix) # strike by p val
      grid.echo()
      plot2 <- grid.grab() %>% as.grob
      matrix.colors <- getGrob(plot2, gPath("circle"), grep = TRUE)[["gp"]][["fill"]]
      # a few adjustments https://stackoverflow.com/questions/53734543/converting-corrplot-output-to-grob
      # lista: childNames(plot2)
      plot2 <- editGrob(plot2,gPath = gPath("circle"), grep = TRUE,
                        gp = gpar(col = NA,
                                  fill = NA))
      plot2 <- editGrob(plot2,
                     gPath("symbols-rect-1"), grep = TRUE,
                     gp = gpar(fill = matrix.colors))

      plot2 <- editGrob(plot2,
                     gPath("background"), grep = TRUE,
                     gp = gpar(fill = NA))
      
      
      # este son las cruces del p value
      plot2 <- editGrob(plot2,
                        gPath("graphics-plot-1-points-1"), grep = TRUE,
                        gp = gpar(fill = NA))
      
      # plot2 <- editGrob(plot2,
      #                   gPath("graphics-plot-1-text-1"), grep = TRUE,
      #                   gp = gpar(col="",fill = NA))
      
      plot2 <- plot2 %>% as.ggplot()
      
      
      ## Convert to graph
      ## ===========================================================================
      # filter by pval
      real_corrs <- real_corrs[real_corrs["pvalue"]<mypval,]
      
      # remove 0 values 
      real_corrs <- real_corrs[real_corrs["weight"]!=0,]
      
      # rename real_corrs (for better plotting afterwards)
      real_corrs[[1]] <- taxa[real_corrs[[1]],] %>% renamer(width=6) %>% gsub(pattern="\\.\\d",replacement = "") 
      real_corrs[[2]] <- taxa[real_corrs[[2]],] %>% renamer(width=6) %>% gsub(pattern="\\.\\d",replacement = "") 

      # create graph
      gr<-igraph::graph_from_data_frame(real_corrs%>%as.matrix, directed = F) 
      
      # color by sign AND weight. (fixed 10 nov para que acepte neg values)
      col1 <- colorRamp(c_scale[6:9])
      col2 <- colorRamp(rev(c_scale[6:9]))
      # E(gr)$colores = apply(col1(E(gr)$weight), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255) )
      v <- c()
      for (w in E(gr)$weight%>%as.numeric()) {
        if (w > 0) {
          x = col1(w)
          v <- c(v, rgb(x[1]/255,x[2]/255,x[3]/255))
        } else {
          x = col2(-w)
          v <- c(v, rgb(x[1]/255,x[2]/255,x[3]/255))
        }
      }
      if (!is.null(v)) { # solo si han sobrevivido nodos al filtro.
        E(gr)$colores <- v
        
        # create layout 
        fr.all <- layout_as_star(gr)
        fr.all <- layout.fruchterman.reingold(gr)
        fr.all.df <- as.data.frame(fr.all)

        # # necesito un pcgcode ordenado
        # pcgcode= pcgcodemaker(taxa[unique(c(real_corrs[[1]],real_corrs[[2]]))], f_ = T)    # TODO maybe entender como lo hice la otra vez
        # mylayouts <- list(GroupByPCG1(as.numeric(colorcode)),
        #             GroupByPCG2(as.numeric(colorcode)))
        # fr.all.df <- as.data.frame(mylayouts[1])
        
    　　# add labels
        fr.all.df$OTUname <- V(gr)$name
        
        # get the edge information using the get.data.frame function
        g <- get.data.frame(gr)
        
        g$from.x <- fr.all.df$V1[match(g$from, fr.all.df$OTUname)]  #  match the from locations from the node data.frame we previously connected
        g$from.y <- fr.all.df$V2[match(g$from, fr.all.df$OTUname)]
        g$to.x <- fr.all.df$V1[match(g$to, fr.all.df$OTUname)]  #  match the to locations from the node data.frame we previously connected
        g$to.y <- fr.all.df$V2[match(g$to, fr.all.df$OTUname)]
        
        # get my graph plotted (ideas from: https://chrischizinski.github.io/rstats/igraph-ggplotll/)
        plot1 <- ggplot() +
          geom_segment(data=g,aes(x=from.x,xend = to.x, y=from.y,yend = to.y,size=weight),colour=g$colores) +
          geom_point(data=fr.all.df,aes(x=V1,y=V2),size=21,colour="black") +  # adds a black border around the nodes
          geom_point(data=fr.all.df,aes(x=V1,y=V2),size=20,colour="lightgrey") +
          geom_text(data=fr.all.df,aes(x=V1,y=V2,label=OTUname)) + # add the node labels
          scale_x_continuous(expand=c(0,250))+  # expand the x limits
          scale_y_continuous(expand=c(0,250))+ # expand the y limits
          theme_void() + # use the ggplot black and white theme
          theme(legend.position = "none",
                panel.background = element_rect(colour = "white", fill = "white")) 
        
      } else {
        # to avoid plotting error                                                # TODO ugly code.
        g <- setNames(data.frame(matrix(ncol=9)), c("from","to","weight","OTUname","from.x","from.y","to.x","to.y", "colores"))
        
        fr.all.df　<- setNames(data.frame(matrix(ncol=3)),c("V1","V2", "OTUname"))
        
        plot1 <- ggplot() +
          geom_segment(data=g,aes(x=from.x,xend = to.x, y=from.y,yend = to.y,size=weight),colour=g$colores) +
          geom_point(data=fr.all.df,aes(x=V1,y=V2),size=21,colour="black") +  # adds a black border around the nodes
          geom_point(data=fr.all.df,aes(x=V1,y=V2),size=20,colour="lightgrey") +
          # geom_text(data=fr.all.df,aes(x=V1,y=V2,label=OTUname)) + # add the node labels
          # scale_x_continuous(expand=c(0,1))+  # expand the x limits 
          # scale_y_continuous(expand=c(0,1))+ # expand the y limits
          theme_void() + # use the ggplot black and white theme
          theme(legend.position = "none",
                panel.background = element_rect(colour = "white", fill = "white"))
      }
      
      # Plot the graph by itself
      ggsave(filename=paste0(outputdir,"/graph_",m,".png"), 
             plot1 + 
               facet_grid(. ~ paste0("Correlation network (",m,")")) +
               theme(strip.background = element_rect(colour="white",fill="white"),
                     strip.text = element_text(size=15, colour="black")),
             height = 7, width = 7)
      
      # Los dos juntos
      # ==============
      bothgg <- ggarrange(
        plot1, NULL,  plot2,
        nrow = 1, widths = c(1, -0.2, 1),
        labels = c(paste0("Correlation network (",m,")"), "", "corrplot with p-values")
      )
      ggsave(bothgg, bg = "white",filename = paste0(outputdir,"/corrs_",m,".png"),width = 15, height=9)
    },
    
    # Si no hay datos de esa m_inic
    error=function(cond){
      message(paste0("Algo fue mal con las muestras del inóculo ",m,", seguramente no hay archivo de simulaciones porque faltaban datos. Mensaje de error:"))
      message(cond)
    }
  )
}
  
  

# Esta es la disrtibución. Se parece mucho a las que salen en el doc suplementario.
# simul_corrs[c,3:n]%>%as.numeric()%>%sort%>%hist()    
# ---> percentil a partir del ecdf
# ---> mirar mis apuntes de teoria de grafos?
# ---> mirar su supplementary
# El problema es que mis p valores son 0 siempre y bueno es porque estoy chequeando una normal. Como chequeo aqui?

# area under the curve to the right of the test statistic z

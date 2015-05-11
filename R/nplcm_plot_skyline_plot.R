#' Distribution of individual prediction for each pathogen sort by specified order
#' 
#' The so-called "skyline-plot" are histograms of individual prediction for each pathogen. 
#' 
#'
#' @param DIR_NPLCM The file path to the folder that stores results from npLCM fit.
#' @param DIR_pathogen_displayorder_lookup The directory path to the .csv file
#'that stores the display order of pathogens in the combined music sheet plot.
#'
#' @return A list of histograms of individual prediction for each pathogen, overall population
#' etiolgy are put in the top-middle section with abbreviated pathogen name; etiology 
#' contributed by each quartile are put at the corresponding section with proportion 
#' contributed to overall etiology by each quartile listed below with brackets. Black histograms
#' are for subset of kids with that particular pathogen detected in NP swab. 
#'
#' @export
#' 
nplcm_plot_skyline_plot <- function(DIR_NPLCM,DIR_pathogen_displayorder_lookup){#BEGIN of function:
            # remember that the data.txt file in winbugs working folder is transposed:
            bugs.dat <- dget(paste(DIR_NPLCM,"data.txt",sep="/"))
            for (bugs.variable.name in names(bugs.dat)) {
              if (!is.null(dim(bugs.dat[[bugs.variable.name]]))) {
                dim(bugs.dat[[bugs.variable.name]]) <- rev(dim(bugs.dat[[bugs.variable.name]]))
                bugs.dat[[bugs.variable.name]] <- aperm(bugs.dat[[bugs.variable.name]])
              }
              assign(bugs.variable.name, bugs.dat[[bugs.variable.name]])
            }
            # change to model_options
            model_options  <- dget(paste(DIR_NPLCM,"model_options.txt",sep="/"))
            pathogen_list <- c(model_options$pathogen_list,model_options$pathogen_SSonly_list)

            VarName ="Icat"

            MBS   = bugs.dat$MBS
            MSS   = bugs.dat$MSS
            MSS.only = bugs.dat$MSS.only
            Jfull_BrS = ncol(MBS)
            JSS_only = ncol(MSS.only)

            Nd = bugs.dat$Nd
            Nu = bugs.dat$Nu
            Y = c(rep(1,Nd),rep(0,Nu))

            MBS.case = bugs.dat$MBS[Y==1,]
            MBS.case.df=as.data.frame(MBS.case)
            SS_index  <- which(colMeans(is.na(bugs.dat$MSS))<.9)
            JSS       <- length(SS_index)
            MSS.case = MSS[,SS_index]
            
            Jfull= Jfull_BrS + JSS_only
            ## reading nplcm outputs:
            res_nplcm <- read.coda(paste(DIR_NPLCM,"coda1.txt",sep="/"),
                                   paste(DIR_NPLCM,"codaIndex.txt",sep="/"),quiet=TRUE)

            # prepare results for individual diagnosis:
            d1_nplcm  = nrow(res_nplcm)
            d2_nplcm  = Nd
            Icat_nplcm = array(NA,c(d1_nplcm,d2_nplcm))
            for (i in 1:d2_nplcm){
              SubVarName = paste(VarName,"[",i,"]",sep="")
              Icat_nplcm[,i] = res_nplcm[,SubVarName]
            }
            # Prepare individual diagnosis table 
            path.prop = apply(Icat_nplcm, 2, function(x) prop.table(table(x)))
            
            freqs.list <- mapply(data.frame,Words=seq_along(path.prop),path.prop,
                                 SIMPLIFY=FALSE,MoreArgs=list(stringsAsFactors=FALSE))
            freqs.df <- do.call(rbind,freqs.list)
            res <- reshape(freqs.df,timevar="Words",idvar="x",direction="wide")
            x.levels = (as.numeric(levels(res$x)))
            row.names(res)=sapply(x.levels,function(x) pathogen_list[x])
            res.sort=res[order(as.numeric(levels(res$x))),]
            res.ggplots = data.frame(t(data.matrix(res.sort)[,-1]))
            
            
            # merge with pathogen_list and make sure every pathogen get a column
            if (ncol(res.ggplots)<length(pathogen_list)){
            path.df = as.data.frame(t(array(0,length(pathogen_list))))
            colnames(path.df)=pathogen_list
            res.ggplots = rbind.fill(path.df,res.ggplots)[-1,]
            }
            
            res.ggplots[is.na(res.ggplots)]=0
            # add MBS to individual diagnosis
            colnames(MBS.case.df)=sapply(colnames(res.ggplots)[1:Jfull_BrS],function(a) paste0(a,"_NPPCR"))
            
            #Find appropriate order
            pathogen_displayorder_lookup <- read.csv(DIR_pathogen_displayorder_lookup)
            f <- pathogen_displayorder_lookup$Pathogen
            display_order <- as.character(levels(f))[f]
            
            pathogen_ord <- rep(NA,length(pathogen_list))
            incre <- 0
            for (j in seq_along(display_order)){
              if (display_order[j]%in% pathogen_list){
                incre <- incre + 1
                pathogen_ord[incre] <- which(pathogen_list==display_order[j])
              }
            }
            #number of each categories
            temp_cat_ind <- sapply(model_options$pathogen_list,
                                   function(path) {which(model_options$pathogen_cat$X==path)})
            temp_cat     <- model_options$pathogen_cat[temp_cat_ind,]
            n.virus = length(temp_cat$pathogen_type[temp_cat$pathogen_type=='V'])
            n.other = 3
            n.bac = length(pathogen_list)-n.virus-n.other
            #re-order the plotting dataset
            res.ggplots = res.ggplots[pathogen_ord]
            #Bacteria vs viruses
            #WF it's not generic
            res.ggplots$VIRUS = rowSums(res.ggplots[,1:n.virus])
            res.ggplots$BACTERIA = rowSums(res.ggplots[,(n.virus+1):(n.bac+n.virus)])
            res.ggplots$OTHER = rowSums(res.ggplots[,(n.bac+n.virus+1):length(pathogen_list)])
            
            #WF export individual probability for each kids
            write.table(res.ggplots, paste0(DIR_NPLCM,"/ind_prop.txt"), sep="\t")
            #second layer of histogram
            res.ggplots_2=cbind(res.ggplots,MBS.case.df)
            #
            res.weight = as.data.frame(array(NA,c(length(pathogen_list)+3,1)))
            res.weight$p4=apply(res.ggplots,2, function(a) paste0(round(100*sum(a[a>=.75])/length(a),2),"%"))
            res.weight$p3=apply(res.ggplots,2, function(a) paste0(round(100*sum(a[a>=.5 & a<.75])/length(a),2),"%"))
            res.weight$p2=apply(res.ggplots,2, function(a) paste0(round(100*sum(a[a>=.25 & a<.5])/length(a),2),"%"))
            res.weight$p1=apply(res.ggplots,2, function(a) paste0(round(100*sum(a[a<.25])/length(a),2),"%"))
            res.weight$w4=apply(res.ggplots,2, function(a) paste0("(",round(100*sum(a[a>=.75])/sum(a),2),"%)"))
            res.weight$w3=apply(res.ggplots,2, function(a) paste0("(",round(100*sum(a[a>=.5 & a<.75])/sum(a),2),"%)"))
            res.weight$w2=apply(res.ggplots,2, function(a) paste0("(",round(100*sum(a[a>=.25 & a<.5])/sum(a),2),"%)"))
            res.weight$w1=apply(res.ggplots,2, function(a) paste0("(",round(100*sum(a[a<.25])/sum(a),2),"%)"))
            res.weight$p=apply(res.ggplots,2, function(a) paste0(round(100*mean(a),2),"%"))
            
            res.weight=res.weight[,-1]
            rownames(res.weight)=colnames(res.ggplots)
          

            #head(res)
            p = list()
            for (i in 1:Jfull){
              name= pathogen_list[pathogen_ord[i]]
              name_NPPCR = paste0(name,"_NPPCR")
              p[[i]]=ggplot(res.ggplots, aes_string(x=name)) + geom_histogram(binwidth = 0.02,colour="black", fill="white") + xlim(0,1.04) +
                theme(axis.text.x  = element_text(angle=0, size=16,face="bold"),
                      axis.text.y  = element_text(angle=90, size=16,face="bold"),
                      axis.title.x  = element_text(size=16,face="bold"),
                      axis.title.y  = element_text(size=20,face="bold"))
              y.max = max(ggplot_build(p[[i]])$data[[1]]$count)
              if (i<=Jfull_BrS-1 | name=="PCP"){
              p[[i]]=p[[i]]+ geom_histogram(data=subset(res.ggplots_2,eval(as.symbol(name_NPPCR))==1),binwidth = 0.02,colour="black", fill="black")
              }
              p[[i]]=p[[i]]+annotate("text", label = res.weight[name,]$w1, x = .125,y = y.max*8/12, size = 9, colour = "black")+
                annotate("text", label = res.weight[name,]$w2, x = .375,y = y.max*8/12, size = 9, colour = "black")+
                annotate("text", label = res.weight[name,]$w3, x = .625,y = y.max*8/12, size = 9, colour = "black")+
                annotate("text", label = res.weight[name,]$w4, x = .875,y = y.max*8/12, size = 9, colour = "black")+
              annotate("text", label = res.weight[name,]$p1, x = .125,y = y.max*9/12, size = 9, colour = "black")+
                annotate("text", label = res.weight[name,]$p2, x = .375,y = y.max*9/12, size = 9, colour = "black")+
                annotate("text", label = res.weight[name,]$p3, x = .625,y = y.max*9/12, size = 9, colour = "black")+
                annotate("text", label = res.weight[name,]$p4, x = .875,y = y.max*9/12, size = 9, colour = "black")+
                annotate("text", label = paste0(name,":",res.weight[name,]$p),x = .5,y = y.max*11/12, size = 13, colour = "black")+
                annotate("rect", xmin = .3, xmax = .7, ymin = y.max*21/24, ymax =y.max* 23/24,alpha = .2)
            }
            p
            
            
}#END of function


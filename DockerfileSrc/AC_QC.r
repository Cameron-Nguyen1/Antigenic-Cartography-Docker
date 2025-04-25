#!/bin/env Rscript
dep = c("ggplot2","stringr","reshape2","RColorBrewer","paletteer","Racmacs","methods","optparse")
lapply(dep,library,character.only=TRUE)

option_list = list(
    make_option(c("--file"), type="character", default=NULL, dest="file",
        help="Titer Table.", metavar="character"),
    make_option(c("--prefix"), type="character", default="MY_EXPERIMENT", dest="prefix",
        help="Meaningful name to apply to output.", metavar="character"),
    make_option(c("--cpu"), type="character", default="2", dest="cpu",
        help="Meaningful name to apply to output.", metavar="character"),
    make_option(c("--repsxdim"), type="character", default="1000", dest="repsxdim",
        help="Number of replicates per dimension, default is 1000.", metavar="character")
)
p = parse_args(OptionParser(option_list=option_list))

options(RacOptimizer.num_cores = parallel::detectCores(as.numeric(p[3])))
set.seed(8675309)
dtest = function(file,prefix,repsxdim){
    y = data.frame(read.csv(as.character(file)))
    nantigen = ncol(y)-1
    nsera = nrow(y)
    rownames(y) = paste0("S_",y[,1])
    y = y[,-1]
    y[1:nantigen] = lapply(y[1:nantigen],as.numeric)
    y = t(y)

    x_log2 = log2(y/10)

    #Caluclate Distance Table using Log Titer Table
    vec = list()
    for (i in 1:ncol(x_log2)){
        vec = append(vec, max(x_log2[,i]))
    }
    for (i in 1:length(vec)){
    x_log2[,i] = vec[[i]]  - x_log2[,i]
    }

    map = acmap(titer_table = y)
    agNames(map) = row.names(x_log2)
    agGroups(map) = row.names(x_log2)
    dilutionStepsize(map) = 2
    map = optimizeMap(map = map, number_of_dimensions=2, number_of_optimizations=1000, minimum_column_basis="none")

    mapdist = data.frame(mapDistances(map))
    tabledist = data.frame(tableDistances(map))
    q = data.frame("Table","Map")
    for (i in 1:nrow(tabledist)){
        for (j in 1:ncol(tabledist)){
        q = rbind(q,c(tabledist[i,j],mapdist[i,j]))
        }
    }
    q=q[-1,]
    colnames(q) = c("Dist_Table","Dist_Map")
    q = data.frame(lapply(q,as.numeric))
    model = lm(Dist_Map~Dist_Table+0,data=q)
    sum = summary(model)
    rval = cor(q$Dist_Table,q$Dist_Map,method="pearson")
    adjrval=sum$adj.r.squared
    qcp = ggplot(q, aes(x=Dist_Table,y=Dist_Map))+ ##Plot
        ggtitle(paste0(prefix," Regresssion, Map Distance ~ Table Distance"))+
        geom_text(x=((max(q$Dist_Table)/2)+1),y=1.75,label=paste0("y = ",signif(model$coefficients[[1]],4),"x + 0"))+
        geom_text(x=((max(q$Dist_Table)/2)+1),y=1.25,label=paste0("R.Val = ",signif(rval,4)))+
        geom_text(x=((max(q$Dist_Table)/2)+1),y=.75,label=paste0("Adj.R^2.Val = ",signif(adjrval,4)))+
        geom_point()+
        stat_smooth(method = lm,formula=y~x+0) +
        #geom_line(aes(y=model2_prediction$lwr),col="coral2",linetype="dashed")+
        #geom_line(aes(y=model2_prediction$upr),col="coral2",linetype="dashed")+
        ylab("Map Distances")+
        xlab("Table Distances")
    ggsave(filename=paste0(prefix,"_Regression.pdf"),plot=qcp,width=10,height=6)

    Dtest = dimensionTestMap(map,replicates_per_dimension = as.numeric(repsxdim))
    write.csv(Dtest,file=paste0(prefix,"_Dimension_Test_RMSE.csv"))

}

dtest(file=p[1],prefix=p[2],repsxdim=p[4])
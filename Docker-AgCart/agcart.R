library(optparse)
option_list = list(
  make_option(c("--Input"), type="character", default=NULL,help="A CSV containing neutralization data with Antigens as Columns and Sera as Rows.", metavar="character"),
  make_option(c("--XY_Lim"), type="character", default="-10,10,-10,10",help="X and Y limits of Cartography graph. Example: \"-6,6,-6,6\" or \"-10,10,-4,4\"", metavar="character"),
  make_option(c("--Prefix"), type="character", default="MyStudy",help="A prefix prepended to generated file names. \"NCBI_SOM\" would give you \"NCBI_SOM_Distance_Between_Groups\".", metavar="character"),
  make_option(c("--Out"), type="character", default="default_output",help="Directs output to designated directory.", metavar="character"),
  make_option(c("--Point_Sizes"), type="character", default="5,2",help="Point sizes of Antigen and Sera coordinates on the AC Map. Antigen,Sera. Example: \"5,2\" or \"4.5,1.5\"", metavar="character"),
  make_option(c("--Transparency"), type="character", default=".8,.4",help="Point transparency of Antigen and Sera coordinates on the AC Map. Antigen,Sera. Example: \".9,.5\" or \".7,.4\"", metavar="character"),
  make_option(c("--Antigen_Overprint"), type="character", default="FALSE",help="Paint Antigens over Sera? TRUE or FALSE.", metavar="character"),
)
p = parse_args(OptionParser(option_list=option_list))
dep = c("ggplot2","stringr","reshape2","RColorBrewer","paletteer","Racmacs")
lapply(dep,library,character.only=TRUE)

dir.create(as.character(p[4]))
y = data.frame(read.csv(file=as.character(p[1])))
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
a = paste0(getwd(),"/",p[4],"/",p[3],"_Cart_Distance_between_groups.csv")
a = as.character(a)
write.csv(x_log2, a)

#Create Antigenic Map
map = acmap(titer_table = y)
agNames(map) = row.names(x_log2)
agGroups(map) = row.names(x_log2)
dilutionStepsize(map) = 2
set.seed(123)
map = optimizeMap(map = map, number_of_dimensions=2, number_of_optimizations=1000, minimum_column_basis="none")

#Save coords to calculate distance of Antigen-Antigen and Sera-Sera
save.coords(map, paste0(p[4],"/",p[3],"_Cart_racmacs_coordinates.csv"))

#Read coord data to calculate within-group distances!
frame68 = read.csv(paste0(p[4],"/",p[3],"_Cart_racmacs_coordinates.csv"))
df2 = frame68[str_detect(frame68$type, "antigen"),] #ANTIGEN DISTANCES
df2 = subset(df2, select = -c(1))
rownames(df2) = df2[,1]
df2 = subset(df2, select = -c(1))
AgDistances = reshape2::melt(as.matrix(dist(df2))) #Distance Matrix command
colnames(AgDistances) = c("Virus 1","Virus 2","Distance")
write.csv(AgDistances,paste0(p[4],"/",p[3],"_Cart_Distance_within_ag.csv"))

df3 = frame68[str_detect(frame68$type, "sera"),] #SERA DISTANCES
rownames(df3) = df3$name
df3 = df3[,c(-1,-2)]
SeraDistances = reshape2::melt(as.matrix(dist(df3)))
colnames(SeraDistances)=c("Sera 1", "Sera 2", "Distance")
write.csv(SeraDistances,paste0(p[4],"/",p[3],"_Cart_Distance_within_sera.csv"))

alpha = str_split_fixed(p[6], ",", 2) #Alphas

#Plots'n'stuff
frame68$color = c(paletteer_d("ggthemes::calc")[1:nantigen],rep("#0d0101",nsera)) #Supports up to 12 Distinct Colors for Antigens, approximately black color applied to Sera.
frame68$size= c(rep(as.numeric(1.5),nantigen),rep(as.numeric(4.5),nsera))
frame68$type = c(rep("Antigen",nantigen),rep("Sera",nsera))
frame68$name = c(rep("Sera",nantigen),frame68$name[112:122])
frame68$alpha = c(rep(as.numeric(.4),nantigen),rep(as.numeric(.8),nsera))
z = "-10,10,-10,10" #XY_lims
#if (toupper(p[7]) == "TRUE"){
#    frame68 = rbind(frame68[grepl("Sera",frame68$type),],frame68[grepl("Antigen",frame68$type),])
#}

pdf(paste0(p[4],"/",p[3],"_",z[2],"x",z[4],".pdf"))
mp =ggplot(frame68, aes(x=X,y=X.1)) +
    geom_point(alpha=frame68$alpha,size=frame68$size,
               aes(shape = type, color = name)) +
    scale_color_manual("Name", values=c(frame68[grepl("Antigen",frame68$type),5],frame68$color[1]),breaks=frame68[grepl("Antigen",frame68$type),2])+
    ylab("Y") +
    coord_fixed(ylim=c(as.numeric(z[1]),as.numeric(z[2])),xlim=c(as.numeric(z[3]),as.numeric(z[4])))+
    scale_x_continuous(breaks=scales::pretty_breaks(n=(as.numeric(z[2])*2)+1))+
    scale_y_continuous(breaks=scales::pretty_breaks(n=(as.numeric(z[4])*2)+1))+
    theme(panel.grid.minor= element_blank())
mp
dev.off()
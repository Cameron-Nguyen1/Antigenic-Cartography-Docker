
set.seed(8675309)
dep = c("ggplot2","stringr","reshape2","paletteer","Racmacs","optparse")
lapply(dep,library,character.only=TRUE)
option_list = list(
  make_option(c("--input"), type="character", default=NULL,help="A CSV containing neutralization data with Antigens as Columns and Sera as Rows.", metavar="character"),
  make_option(c("--xy_lim"), type="character", default="-10,10,-10,10",help="X and Y limits of Cartography graph. Example: \"-6,6,-6,6\" or \"-10,10,-4,4\"", metavar="character"),
  make_option(c("--prefix"), type="character", default="MyStudy",help="A prefix prepended to generated file names. \"NCBI_SOM\" would give you \"NCBI_SOM_Distance_Between_Groups\".", metavar="character"),
  make_option(c("--psizes"), type="character", default="5,2",help="Point sizes of Antigen and Sera coordinates on the AC Map. Antigen,Sera. Example: \"5,2\" or \"4.5,1.5\"", metavar="character"),
  make_option(c("--opacity"), type="character", default=".8,1",help="Point transparency of Antigen and Sera coordinates on the AC Map. Antigen,Sera. Example: \".9,.5\" or \".7,.4\"", metavar="character"),
  make_option(c("--agoverprint"), type="character", default="FALSE",help="Paint Antigens over Sera? TRUE or FALSE.", metavar="character"),
  make_option(c("--agsort"), type="character", default="TRUE",help="Sort Antigens? TRUE or FALSE.", metavar="character")
)

pargs = parse_args(OptionParser(option_list=option_list))

y = data.frame(read.csv(file=as.character(pargs$input)))
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
a = paste0(pargs$prefix,"_Cart_Distance_between_groups.csv")
a = as.character(a)
write.csv(x_log2, a)

#Create Antigenic Map
map = acmap(titer_table = y)
agNames(map) = row.names(x_log2)
agGroups(map) = row.names(x_log2)
dilutionStepsize(map) = 2
map = optimizeMap(map = map, number_of_dimensions=2, number_of_optimizations=1000, minimum_column_basis="none")
saveRDS(map,paste0(pargs$prefix,"_mapobj.rds"))

#Save coords to calculate distance of Antigen-Antigen and Sera-Sera
save.coords(map, paste0(pargs$prefix,"_Cart_racmacs_coordinates.csv"))

#Read coord data to calculate within-group distances!
df1 = read.csv(paste0(pargs$prefix,"_Cart_racmacs_coordinates.csv"))
if (toupper(pargs$agsort) == "TRUE"){
  df1[1:nantigen,] = df1[1:nantigen,][order(df1[1:nantigen,]$name),]
}
df2 = df1[str_detect(df1$type, "antigen"),] #ANTIGEN DISTANCES
df2 = subset(df2, select = -c(1))
rownames(df2) = df2[,1]
df2 = subset(df2, select = -c(1))
AgDistances = reshape2::melt(as.matrix(dist(df2))) #Distance Matrix command
colnames(AgDistances) = c("Virus 1","Virus 2","Distance")
write.csv(AgDistances,paste0(pargs$prefix,"_Cart_Distance_within_ag.csv"))

df3 = df1[str_detect(df1$type, "sera"),] #SERA DISTANCES
rownames(df3) = df3$name
df3 = df3[,c(-1,-2)]
SeraDistances = reshape2::melt(as.matrix(dist(df3)))
colnames(SeraDistances)=c("Sera 1", "Sera 2", "Distance")
write.csv(SeraDistances,paste0(pargs$prefix,"_Cart_Distance_within_sera.csv"))

alpha = str_split_fixed(pargs$opacity, ",", 2) #Alphas

#Plots'n'stuff
if (nantigen > 20)
  {
    df1$colors = c(paletteer_d("ggsci::default_igv")[1:nantigen],rep("#0d0101",nsera))
  }else{
    df1$colors = c(paletteer_d("ggsci::category20_d3")[1:nantigen],rep("#0d0101",nsera))
}
size = str_split_fixed(pargs$psizes,",",2) #Sizes
df1$size= c(rep(as.numeric(size[1]),nantigen),rep(as.numeric(size[2]),nsera))
df1$Type = c(rep("Antigen",nantigen),rep("Sera",nsera))
df1$name = c(df1$name[1:nantigen],rep("Sera",nsera))
df1$alpha = c(rep(as.numeric(alpha[1]),nantigen),rep(as.numeric(alpha[2]),nsera))
z = str_split_fixed(pargs$xy_lim, ",", 4) #XY_lims
if (toupper(pargs$agoverprint) == "TRUE"){
  df1 = rbind(df1[grepl("Sera",df1$Type),],df1[grepl("Antigen",df1$Type),])
}

mypdf = paste0(pargs$prefix,"_Cartography.pdf")
mp =ggplot(df1, aes(x=X,y=X.1)) +
    geom_point(alpha=df1$alpha,size=df1$size,
              aes(shape = Type, color = name)) +
    # form_blobs_2d(blobject=blobject.95,meta=df1) +
    scale_color_manual("Name", values=c(df1[grepl("Antigen",df1$Type),5],"#0d0101"),breaks=df1[grepl("Antigen",df1$Type),2])+
    ylab("Y") +
    coord_fixed(ylim=c(as.numeric(z[1]),as.numeric(z[2])),xlim=c(as.numeric(z[3]),as.numeric(z[4])))+
    scale_x_continuous(breaks=scales::pretty_breaks(n=(as.numeric(z[2])*2)+1))+
    scale_y_continuous(breaks=scales::pretty_breaks(n=(as.numeric(z[4])*2)+1))+
    theme(panel.grid.minor= element_blank(),axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.ticks = element_line(size=1.15),
          legend.text = element_text(size = 12), legend.title = element_text(size = 14))
ggsave(filename=mypdf,plot=mp)

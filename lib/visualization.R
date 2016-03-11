########################################################
## Visualization: MDS, HeatMaps, BoxPlot              ##
########################################################

flog.info("make MDS for different taxonomy levels")
for (i in names(totalTable)[1:3]){
  dist<-bcdist(totalTable[[i]][which(totalTable$meta[,"Type.1"] %in% "EKOK_1"),])
  MDS(dist, pathway$group.one$outdir, totalTable$meta[which(totalTable$meta[,"Type.1"] %in% "EKOK_1"),],
      "Type.1", i)
}

for (i in names(totalTable)[1:3]){
  dist<-bcdist(totalTable[[i]][which(totalTable$meta[,"Type.1"] %in% "EKOK_2"),])
  MDS(dist, pathway$group.one$outdir, totalTable$meta[which(totalTable$meta[,"Type.1"] %in% "EKOK_2"),],
      "Type.1", i)
}

flog.info("make heatmaps")
#Making HeatMap
HeatMap(fam.gen.spe.list.g1$family, pathway$group.one$outdir, 'family.avva.g1' )
HeatMap(fam.gen.spe.list.g1$genus, pathway$group.one$outdir, 'genus.avva.g1')
HeatMap(fam.gen.spe.list.g1$species, pathway$group.one$outdir, 'species.avva.g1')

HeatMap(fam.gen.spe.list.g2$family, pathway$group.one$outdir, 'family.avva.g2')
HeatMap(fam.gen.spe.list.g2$genus, pathway$group.one$outdir, 'genus.avva.g2')
HeatMap(fam.gen.spe.list.g2$species, pathway$group.one$outdir, 'species.avva.g2')

#Making HeatMap for TOP Features
HeatMap(fam.gen.spe.listTOP.g1$family, pathway$group.one$outdir, 'family.avvaTOPfeatures.g1' )
HeatMap(fam.gen.spe.listTOP.g1$genus, pathway$group.one$outdir, 'genus.avvaTOPfeatures.g1' )
HeatMap(fam.gen.spe.listTOP.g1$species, pathway$group.one$outdir, 'species.avvaTOPfeatures.g1')

HeatMap(fam.gen.spe.listTOP.g2$family, pathway$group.one$outdir, 'family.avvaTOPfeatures.g2' )
HeatMap(fam.gen.spe.listTOP.g2$genus, pathway$group.one$outdir, 'genus.avvaTOPfeatures.g2' )
HeatMap(fam.gen.spe.listTOP.g2$species, pathway$group.one$outdir, 'species.avvaTOPfeatures.g2')


  
  
  


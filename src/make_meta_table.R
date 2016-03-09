########################################
# making meta table from otu files     #
########################################
current_directory <- "~/home/metagenome/AVVA"
path.to.meta <- "/home/anna/metagenome/AVVA/src/pathway_meta.ini"
meta.pathway <- ReadIni(path.to.meta) 

Load <- function(pathway)
{
  group.one <- ReadQiimeSumFeats(meta.pathway$group.one$fam.group.one.otu.tbl)
  group.two <- ReadQiimeSumFeats(meta.pathway$group.two$fam.group.tow.otu.tbl)
  group.three <- ReadQiimeSumFeats(meta.pathway$group.three$fam.group.three.otu.tbl)
  group.four <- ReadQiimeSumFeats(meta.pathway$group.four$fam.group.four.otu.tbl)
  group.five <- ReadQiimeSumFeats(meta.pathway$group.five$fam.group.five.otu.tbl)
  group.six <- ReadQiimeSumFeats(meta.pathway$group.six$fam.group.six.otu.tbl)
  list(group.one=group.one, group.two=group.two, group.three=group.three, 
       group.four=group.four, group.five=group.five, group.six=group.six)
}

meta_name <- Load(meta.pathway)

g.one <- as.data.frame(meta_name$group.one)
g.two <- as.data.frame(meta_name$group.two)
g.three <- as.data.frame(meta_name$group.three)
g.four <- as.data.frame(meta_name$group.four)
g.five <- as.data.frame(meta_name$group.five)
g.six <- as.data.frame(meta_name$group.six)


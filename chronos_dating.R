library(ape)

setwd("/home/mlb/Phylo-data/UCEs/Chronos/")

bin1_tree <- read.tree("RAxML_bestTree.ML-bin1-part")

rb1tr <- root(bin1_tree, "Harpegnathos_saltator_Genome", resolve.root=T)

mrca(rb1tr)["Harpegnathos_saltator_Genome", "Acanthostichus_serratulus_M166"]
mrca(rb1tr)["Camponotus_floridanus_Genome", "Atta_cephalotes_Genome"]
mrca(rb1tr)["Pogonomyrmex_barbatus_Genome", "Atta_cephalotes_Genome"]
mrca(rb1tr)["Chrysapace_TH_M156", "Yunodorylus_TH_M191"]
mrca(rb1tr)["Cylindromyrmex_meinerti_D778", "Cylindromyrmex_darlingtoni_M211"]
mrca(rb1tr)["Eciton_hamatum_M293", "Neivamyrmex_californicus_M272"]

calib <- makeChronosCalib(rb1tr, node=165, age.min=94.3, age.max=135, 
                          interactive=F)
calib <- rbind(calib, makeChronosCalib(rb1tr, node=322, age.min=89.3, age.max=120, 
                                       interactive=F))
calib <- rbind(calib, makeChronosCalib(rb1tr, node=323, age.min=33.9, age.max=110, 
                                       interactive=F))
calib <- rbind(calib, makeChronosCalib(rb1tr, node=205, age.min=33.9, age.max=85, 
                                       interactive=F))
calib <- rbind(calib, makeChronosCalib(rb1tr, node=222, age.min=13.7, age.max=43, 
                                       interactive=F))
calib <- rbind(calib, makeChronosCalib(rb1tr, node=235, age.min=16, age.max=46, 
                                       interactive=F))

c <- chronos(rb1tr, lambda = 1, model = "discrete", quiet=F, 
             calibration = calib,
             control = chronos.control(nb.rate.cat=5))
attr(c, "PHIIC")[3]
write.tree(c, file="calib_bin1_outgr.discrete5")
c <- chronos(rb1tr, lambda = 1, model = "discrete", quiet=F, 
             calibration = calib,
             control = chronos.control(nb.rate.cat=1))
attr(c, "PHIIC")[3]
write.tree(c, file="calib_bin1_outgr.strict")
c <- chronos(rb1tr, lambda = 1, model = "correlated", quiet=F, 
             calibration = calib,
             control = chronos.control())
write.tree(c, file="calib_bin1_outgr.correlated1.0")
attr(c, "PHIIC")[4]
c <- chronos(rb1tr, lambda = 0.5, model = "correlated", quiet=F, 
             calibration = calib,
             control = chronos.control())
write.tree(c, file="calib_bin1_outgr.correlated0.5")
attr(c, "PHIIC")[4]
c <- chronos(rb1tr, lambda = 0.0, model = "correlated", quiet=F, 
             calibration = calib,
             control = chronos.control())
attr(c, "PHIIC")[4]
write.tree(c, file="calib_bin1_outgr.correlated0.0")

c <- chronos(rb1tr, lambda = 1, model = "relaxed", quiet=F, 
             calibration = calib,
             control = chronos.control())
attr(c, "PHIIC")[4]
write.tree(c, file="calib_bin1_outgr.relaxed1.0")
c <- chronos(rb1tr, lambda = 0.5, model = "relaxed", quiet=F, 
             calibration = calib,
             control = chronos.control())
attr(c, "PHIIC")[4]
write.tree(c, file="calib_bin1_outgr.relaxed0.5")
c <- chronos(rb1tr, lambda = 0.0, model = "relaxed", quiet=F, 
             calibration = calib,
             control = chronos.control())
attr(c, "PHIIC")[4]
write.tree(c, file="calib_bin1_outgr.relaxed0.0")

for (i in 1:100) {
  c <- chronos(rb1tr, lambda = 1, model = "discrete", quiet=F, 
               calibration = calib,
               control = chronos.control(nb.rate.cat=1))
  print(attr(c, "PHIIC")[3])
  write.tree(c, file=paste("calib_bin1_outgr.strict", i, sep="-"))
}
for (i in 70:100) {
  c <- chronos(rb1tr, lambda = 1, model = "discrete", quiet=F, 
               calibration = calib,
               control = chronos.control(nb.rate.cat=10))
  print(attr(c, "PHIIC")[3])
  write.tree(c, file=paste("calib_bin1_outgr.10discrete", i, sep="-"))
}

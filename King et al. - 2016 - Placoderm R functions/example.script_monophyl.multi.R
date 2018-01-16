# Will return the probability that a taxon group, or multiple taxon groups, are monophyletic in the posterior sample of trees
# Includes the option to assess monophyly with the exclusion of rogue taxa
# Returns the posterior probability of all the clades being monophyletic, as well as all the trees in which they not, and the generations in which they are not, to allow inspection

setwd("C:/Workspace/king0393/Placoderms/R_monophyly_statistic")
library(caper)
source("monophyl.multi2.R")


#Load the trees
read.nexus("placoderms.trees") -> trees


#define taxon groups in vectors
placoderms <- c("Jagorina_pandora", "Gemuendina_stuertzi", "Austroptyctodus_gardineri", "Bothriolepis_canadensis", "Bothriolepis_sp_Gogo", "Brindabellaspis_stensioi", "Buchanosteus_confertituberculatus", "Campbellodus_decipiens", "Coccosteus_cuspidatus", "Compagopiscis_croucheri", "Cowralepis_mclachlani", "Diandongpetalichthys_liaojiaoshanensis", "Dicksonosteus_arcticus", "Eastmanosteus_calliaspis", "Eurycaraspis_incilis", "Gavinaspis_convergens", "Groenlandaspis_sp_Mt_Howitt", "Holonema_westolli", "Incisoscutum_ritchiei", "Kujdanowiaspis_podolica", "Lunaspis_broili", "Macropetalichthys_rapheidolabis", "Materpiscis_attenboroughi", "Microbrachius_dicki", "Parabuchanosteus_murrumbidgeensis", "Parayunnanolepis_xitunensis", "Pterichthyodes_milleri", "Quasipetalichthys_haikouensis", "Remigolepis_walkeri", "Rhamphodopsis_threiplandi", "Romundina_stellina", "Sigaspis_lepidophora", "Wuttagoonaspis_fletcheri", "Yunnanolepis_sp")
osteichthyans <- c("Dialipina_salgueiroensis", "Ligulalepis_toombsi", "Cheirolepis_canadensis", "Cheirolepis_trailli", "Howqualepis_rostridens", "Mimipiscis_toombsi", "Moythomasia_durgaringa", "Kentuckia_deani", "Osorioichthys_marginis", "Meemannia_eos", "Guiyu_oneiros", "Psarolepis_romeri", "Achoania_jarvikii", "Onychodus_jandemarrai", "Miguashaia_bureaui", "Styloichthys_changae", "Diabolepis_speratus", "Youngolepis_praecursor", "Powichthys_thorsteinssoni", "Porolepis_sp", "Glyptolepis_groenlandica", "Kenichthys_campbelli", "Osteolepis_macrolepidota", "Gogonasus_andrewsae", "Eusthenopteron_foordi")
AC <- c("Acanthodes_bronni", "Brachyacanthus_scutiger", "Brochoadmones_milesi", "Cassidiceps_vermiculatus", "Cheiracanthus_sp", "Climatius_reticulatus", "Culmacanthus_stewarti", "Euthacanthus_macnicoli", "Gladiobranchus_probaton", "Homalacanthus_concinnus", "Ischnacanthus_gracilis", "Kathemacanthus_rosulentus", "Latviacanthus_ventspilsensis", "Lupopsyrus_pygmaeus", "Mesacanthus_mitchelli", "Obtusacanthus_corroconis", "Parexus_recurvus", "Poracanthodes_menneri", "Promesacanthus_eppleri", "Ptomacanthus_anglicus", "Diplacanthus_striatus", "Tetanopsyrus_lindoei_breviacanthias", "Vernicomacanthus_waynensis", "Cladodoides_wildungensis", "Akmonistion_zangerli", "Cobelodus_braincase", "Cladoselache_kepleri_fyleri", "Chondrenchelys_problematica", "Helodus_simplex", "Debeerius_ellefseni", "Doliodus_problematicus", "Hamiltonichthys_mapesi", "Onychoselache_traquari", "Orthacanthus_sp", "Pucapampella_rodrigae", "Tamiobatis_vetustus", "Tristychius_arcuatus")
Janusiscus_osteichthyans <- c(osteichthyans, "Janusiscus_schulzei")
Janusiscus_AC <- c(AC, "Janusiscus_schulzei")
Entelognathus_Janusiscus_osteichthyans <- c(Janusiscus_osteichthyans, "Entelognathus_primordialis")
placoderms_Entelognathus_Janusiscus_osteichthyans <- c(Entelognathus_Janusiscus_osteichthyans, placoderms)
Entelognathus_placoderms <- c(placoderms, "Entelognathus_primordialis")
Janusiscus_Entelognathus_placoderms <- c(Entelognathus_placoderms, "Janusiscus_schulzei")
osteichthyans_AC <- c(osteichthyans, AC)
Janusiscus_osteichthyans_AC <- c(osteichthyans_AC, "Janusiscus_schulzei")
Entelognathus_Janusiscus_osteichthyans_AC <- c(Janusiscus_osteichthyans_AC, "Entelognathus_primordialis")
placoderms_AC <- c(placoderms, AC)
Entelognathus_placoderms_AC <- c(placoderms_AC, "Entelognathus_primordialis")
Janusiscus_Entelognathus_placoderms_AC <- c(Entelognathus_placoderms_AC, "Janusiscus_schulzei")


# insert taxon group vectors as separate elements of a list
#Any number including one group is possible

list() -> clades
placoderms -> clades [[1]]
osteichthyans -> clades[[2]]
AC -> clades[[3]]
Janusiscus_osteichthyans -> clades[[4]]
Janusiscus_osteichthyans_AC -> clades[[5]]
Entelognathus_Janusiscus_osteichthyans_AC -> clades[[6]]

# define rogue taxa as a last (to be pruned from all trees within the function)
"Ramirosuarezia_boliviana" -> rogue.taxa

# run it
# The first two options are the trees and the list of clades
# exclude.rogue is whether or not rogue taxa will be excluded. Default false
# if exclude.rogue is TRUE, then a vector of rogue taxa should be included. Rogue taxa should not appear in any of the taxon groups being assessed
# burnin is the burn in. Default 0.1

monophyl.multi(trees, clades, exclude.rogue=T, rogue.taxa, burnin=0.1) -> result
result

# result$posterior is the posterior probability all the clades are monophyltic at the same time
# result$failedtrees are all the trees in which the monophyly conditons are not met. These can be plotted
# result$failedgens are the generations (starting from the beginning, not after burnin) of these failed trees. Useful for looking at the trees in figtree instead



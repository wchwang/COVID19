

library(gprofiler2)
library(magrittr)
library(dplyr)
set_base_url("https://biit.cs.ut.ee/gprofiler_archive3/e98_eg45_p14")
custom_id <- upload_GMT_file("../Data/Reactome/ReactomePathways.gmt")

candidate_drugs<-read.csv("COVID19_candidate_approved_drugs_extended_targets.csv")
#candidate_drugs2<-candidate_drugs %>% filter(name =="Masoprocol" | name =="Unoprostone")

for (row in 1:nrow(candidate_drugs)){
  drug_name<-as.character(candidate_drugs[row,"Name"])
  drug_targets<-as.character(candidate_drugs[row,"Targets"])
  l_drug_targets <- unlist(strsplit(drug_targets,","))
  #drug_targets_chr<-drug_targets_tib$Targets
  print(drug_name)
  print(l_drug_targets)
  print(length(l_drug_targets))
  len_drug_target <- length(l_drug_targets)
  drug_target_enriched_path<-gost(l_drug_targets,organism = custom_id, evcodes = TRUE) # threshold = 0.05
  enriched_result<-drug_target_enriched_path$result
  print(enriched_result)
  drops<-c("parents","evidence_codes")
  enriched_result_df<-enriched_result[,!(names(enriched_result) %in% drops)]
  enrich_result_addr <- sprintf("../Result/Drug_pathways/enriched_reactome_pathways_%s.csv",drug_name)
  print(enrich_result_addr)
  write.csv(enriched_result_df,file = enrich_result_addr)
}


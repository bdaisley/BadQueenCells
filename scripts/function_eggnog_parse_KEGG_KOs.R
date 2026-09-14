
library(dplyr)
library(Biostrings)
library(stringr)
library(reshape2)
library(ggplot2)
library(heatmap3)
library(ggheatmap)
library(readr)
library(ggdendro)
library(cowplot)
library(dendsort)
library(ggiraph)
library(ggordiplots)
library(Maaslin2)
library(ggrepel)

kegg.names <- readr::read_tsv("/home/brendan/.conda/envs/metacerberus/lib/python3.12/site-packages/meta_cerberus/DB/KEGG.tsv", show_col_types=FALSE, col_names=TRUE) %>%
  distinct(ID, .keep_all=TRUE)
#---------------- 1) Read final eggnog mapper annotation files
input.list.egg <- list.files("data/eggnog_results/", pattern=".annotations", 
                             recursive=TRUE)[1:12]
emapper.collector <- data.frame("combined" = "")
for(i in input.list.egg){
  short.name <- str_split_fixed(i, "/", 2)[,1]
  print(short.name)
  emapper.in <- readr::read_tsv(paste0("data/eggnog_results/", i),
                                skip=4, show_col_types=FALSE, col_names=TRUE) %>%
    filter(!grepl("Metazoa|Eukary", max_annot_lvl)) %>% filter(!is.na(seed_ortholog)) %>%
    filter(KEGG_ko!="-") %>%
    mutate(combined=paste0(max_annot_lvl, "___", Description, "___", KEGG_ko)) %>% 
    group_by(combined) %>%
    mutate(counts=n()) %>% dplyr::slice(1) %>%
    dplyr::select(combined, counts) %>% 
    dplyr::rename(!!short.name := counts)
  emapper.collector <- merge(emapper.collector,emapper.in, by="combined", sort=FALSE, all=TRUE)
}

emapper.collector.x <- emapper.collector %>% mutate_at(vars(matches("BCMy")), ~ifelse(is.na(.), 0, .)) %>%
  filter(!grepl("Euk", combined)) %>%
  mutate(KO=str_split_fixed(combined, "___ko:", 2)[,2]) %>%
  dplyr::relocate(KO, .after="combined") %>% dplyr::select(-combined) %>% mutate(KO = gsub("ko:", "", KO)) %>%
  tidyr::separate_rows(KO, sep = ",") %>%
  filter(KO!="") %>%
  group_by(KO) %>% mutate_at(vars(matches("BCMy")), ~sum(.)) %>% dplyr::slice(1) %>% ungroup()

write_tsv(emapper.collector.x, file="data/eggnog_results/eggnog_annotations_parsed.tsv")

#---------------------------------------
# Filter host-associated KOs
#---------------------------------------
eggnog.functions <- read_tsv("data/eggnog_results/eggnog_annotations_parsed.tsv")
emapper.collector.x.mer <- merge(kegg.names, emapper.collector.x, by.x="ID", by.y="KO", sort=FALSE, all=FALSE) %>%
  filter(!grepl("Aging|Cancer|Cardiovascular|Cellular community|Circulatory|Development|Digestive|Endocrine|Excretory|Immune|Infectious disease|Nervous|Neurodegenerative", L2)) %>%
  filter(!grepl("Poorly characterized|Protein families|Sensory|Substance dependence", L2)) %>%
  dplyr::select(matches("ID|BCMy")) %>% dplyr::rename(KO="ID") 

#---------------------------------------
#Make list of KOs for module calculate
#---------------------------------------
for(i in colnames(emapper.collector.x.mer)[2:length(colnames(emapper.collector.x.mer))]){
  emapper.collector.x.mer.filt <- emapper.collector.x.mer %>% dplyr::select(KO, i) %>% filter(get(i) > 0) %>% dplyr::select(KO) %>% sjmisc::rotate_df()
  write_csv(emapper.collector.x.mer.filt, paste0("data/eggnog_results/", i, "/", i, "_KO_list.txt"), col_names = FALSE)
}

#-------------------------------------------------------------
# Run "09_eggnog_05_KEGG_MODULES" script before next step
#-------------------------------------------------------------
collector.mod.df <- data.frame(module_accession="NA")
for(i in colnames(emapper.collector.x.mer)[2:length(colnames(emapper.collector.x.mer))]){
  module.in <- readr::read_tsv(paste0("data/eggnog_results/", i, "/", i, "_KO_list_modules_pathways.tsv"),
                               skip=0, show_col_types=FALSE, col_names=TRUE) %>% dplyr::rename(!!i:= completeness)
  if(i == colnames(emapper.collector.x.mer)[2] ){
    collector.mod.df <- merge(collector.mod.df, module.in, by="module_accession", sort=FALSE, all=TRUE)# %>%
  }
  if (i != colnames(emapper.collector.x.mer)[2] ){
    collector.mod.df <- merge(collector.mod.df, module.in, by="module_accession", sort=FALSE, all=TRUE) #%>%
    #print(colnames(collector.mod.df))
    collector.mod.df <- collector.mod.df %>% mutate(pathway_name.x = ifelse(is.na(pathway_name.x), pathway_name.y, pathway_name.x)) %>%
      mutate(pathway_class.x = ifelse(is.na(pathway_class.x), pathway_class.y, pathway_class.x)) %>%
      mutate(matching_ko.x = ifelse(is.na(matching_ko.x), matching_ko.y, matching_ko.x)) %>%
      mutate(missing_ko.x = ifelse(is.na(missing_ko.x), missing_ko.y, missing_ko.x)) %>%
      dplyr::rename(pathway_name = pathway_name.x, pathway_class = pathway_class.x, matching_ko = matching_ko.x, missing_ko = missing_ko.x) %>%
      dplyr::select(!matches("pathway_name.y|pathway_class.y|matching_ko.y|missing_ko.y"))
    #print(colnames(collector.mod.df))
  }
}
collector.mod.df.final <- cbind(collector.mod.df %>% dplyr::select(!matches("BC")), collector.mod.df %>% dplyr::select(matches("BC"))) %>%
  mutate_at(vars(matches("BC")), ~ifelse(is.na(.), 0, .))

write_tsv(collector.mod.df.final, file="data/eggnog_results/eggnog_annotations_parsed_modules.tsv")


library(ggplot2)
library(broom)
library(tidyverse)
library(tidycat)
source("/projectnb2/talbot-lab-data/zrwerbin/temporal_forecast/functions/helperFunctions.r")



#######################################

# Metagenomic family presence
all_fam_df <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/interactions/fam_abund.rds")

# Families of interest
fam_list <- c("f_Nitratiruptoraceae", "f_Nitrosomonadaceae", "f_Nitrosopumilaceae",
							"f_Nitrososphaeraceae", "f_Nitrospiraceae")
fam_df <- all_fam_df  %>%
	filter(percentage > 0) %>%
	#filter(!plotID %in% c("HARV_010","HARV_021","HARV_033")) %>%
	filter(name %in% fam_list)
fam_df <- fam_df %>% mutate(siteID = substr(sampleID, 1, 4)) %>%
	mutate(dates = sapply(strsplit(sampleID, "-"), function(x) x[grep("[2]\\d\\d\\d\\d\\d\\d\\d", x)])) %>%
	mutate(dateID = substr(as.character(dates), 1, 6)) %>%
	mutate(plotID = substr(sampleID, 1, 8)) %>%
	mutate(site_date = paste0(siteID, "-", dateID)) %>%
	mutate(horizon = ifelse(grepl("-M-", sampleID), "M", "O")) %>%
	mutate(plot_date = paste0(plotID, "-", dateID)) %>%
	as.data.frame()
fam_df_scaled <- fam_df %>% group_by(name) %>% mutate(scaled_percent = scale(percentage))


#######################################

# METAGENOMIC - GENERA OF INTEREST

abun_df <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/interactions/genus_abund.rds")

genera_list <- c("g_Bacillus", "g_Methanosarcina", "g_Pseudomonas",
								 "g_Saccharomyces", "g_Staphylococcus", "g_Clostridium","g_Nitrobacter","g_Nitrospira","g_Nitrospirillum","g_Nitratireductor","g_Geobacter","g_Acinetobacter","g_Klebsiella", "g_Mycobacterium","g_Vibrio","g_Nitrosomonas")
genera_df <- abun_df %>% filter(name %in% genera_list)
genera_df <- genera_df %>% mutate(siteID = substr(sampleID, 1, 4)) %>%
	mutate(dates = substr(sampleID, 20, 27)) %>%
	mutate(asDate = as.Date(as.character(dates), "%Y%m%d")) %>%
	mutate(dateID = substr(as.character(dates), 1, 6)) %>%
	mutate(plotID = substr(sampleID, 1, 8)) %>%
	mutate(site_date = paste0(siteID, "-", dateID)) %>%
	mutate(horizon = ifelse(grepl("-M-", sampleID), "M", "O")) %>%
	mutate(plot_date = paste0(plotID, "-", dateID)) %>%
	as.data.frame()
genera_df_scaled <- genera_df %>% group_by(name) %>% mutate(scaled_percent = scale(percentage))
meta_gen = genera_df_scaled %>% filter(name %in% c("g_Nitrosomonas","g_Nitrospira","g_Pseudomonas","g_Geobacter"))


#######################################


# 16S - GENERA OF INTEREST
nitr_abun2 <- readRDS("/projectnb2/talbot-lab-data/zrwerbin/interactions/data/nitrifiers_genus_16S.rds")
nitr_abun2 <- nitr_abun2 %>% mutate(nitrifier_genera = candidatus.nitrosotalea + candidatus.nitrotoga + nitrolancea +
																			nitrosomonadaceae_genus + nitrosospira + nitrospira) #%>% select(nitrifier_genera)
#nitr_abun2$geneticSampleID <- rownames(nitr_abun2)
nitr_abun2$sampleID = gsub("-DNA[123]","",rownames(nitr_abun2))

nitr_abun2 <- nitr_abun2 %>% mutate(siteID = substr(sampleID, 1, 4)) %>%
	mutate(geneticSampleID = sapply(strsplit(sampleID, "-DNA"),  "[[" , 1)) %>%
	mutate(sampleID = sapply(strsplit(sampleID, "-gen.fastq"),  "[[" , 1)) %>%
	mutate(dates = sapply(strsplit(sampleID, "-"), function(x) x[grep("[2]\\d\\d\\d\\d\\d\\d\\d", x)])) %>%
	mutate(dates = ifelse(dates == "21040514", "20140514", dates)) %>%
	mutate(asDate = as.Date(as.character(dates), "%Y%m%d")) %>%
	mutate(dateID = substr(as.character(dates), 1, 6)) %>%
	mutate(plotID = substr(sampleID, 1, 8)) %>%
	mutate(site_date = paste0(siteID, "-", dateID)) %>%
	mutate(horizon = ifelse(grepl("-M-", sampleID), "M", "O")) %>%
	mutate(plot_date = paste0(plotID, "-", dateID)) %>%
	as.data.frame()
nitr_gen_abun <- nitr_abun2 %>% filter(siteID=="HARV")



#######################################


# 16S - FUNCTIONAL GROUPS
# Read in abundances from pathway presence ("functional groups")
cal <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds")
val <- readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/val_groupAbundances_16S_2021.rds")

nitr_abun <- rbind(cal$nitrification, val$nitrification)
nitr_abun$geneticSampleID <- gsub("-DNA[123]","",nitr_abun$sampleID)
nitr_fg_abun <- nitr_abun %>%
	mutate(horizon = ifelse(grepl("-M-", sampleID), "M", "O")) %>%
	rename(nitrifier_abundance = nitrification) %>%
	filter(siteID=="HARV")



# OTHER 16S - GENERA OF INTEREST

gen_abun <- rbind(cal$genus_bac, val$genus_bac)

gen_abun$geneticSampleID <- gsub("-DNA[123]","",gen_abun$sampleID)
gen_abun <- gen_abun %>% mutate(horizon = ifelse(grepl("-M-", sampleID), "M", "O")) %>%
	rename(nitrospira_16S_abundance = nitrospira) %>%  filter(siteID=="HARV")
# THIS ONE GOOD
gen_abun %>% filter(siteID=="HARV") %>% ggplot(aes(x = horizon, y = nitrospira_16S_abundance)) + geom_boxplot() + geom_jitter(alpha=.5)


#######################################

# 16S - FAMILIES OF INTEREST

family_abun <- rbind(cal$family_bac, val$family_bac)





# misc???

other_fam_df <- all_fam_df %>%
	#filter(!plotID %in% c("HARV_010","HARV_021","HARV_033")) %>%
	filter(name %in% c("f_Pseudomonadaceae","f_Methanosarcinaceae"))
other_fam_df <- other_fam_df %>% mutate(siteID = substr(sampleID, 1, 4)) %>%
	mutate(dates = sapply(strsplit(sampleID, "-"), function(x) x[grep("[2]\\d\\d\\d\\d\\d\\d\\d", x)])) %>%
	mutate(dateID = substr(as.character(dates), 1, 6)) %>%
	mutate(plotID = substr(sampleID, 1, 8)) %>%
	mutate(site_date = paste0(siteID, "-", dateID)) %>%
	mutate(horizon = ifelse(grepl("-M-", sampleID), "M", "O")) %>%
	mutate(plot_date = paste0(plotID, "-", dateID)) %>%
	as.data.frame()
other_fam_df_scaled <- other_fam_df %>% group_by(name) %>% mutate(scaled_percent = scale(percentage))
other_fam_df_scaled %>%
	ggplot(aes(x = horizon, y = percentage, color = name)) +
	geom_boxplot(position="dodge") +
	geom_jitter(alpha=.5)






# Add sample IDs
eval_df <- read.csv("/projectnb/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/true_rates_df.csv")

biomass_long <- readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/extramods_total_biomass.rds")
biomass_long <- readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/more_sugars_total_biomass.rds")
biomass_long <- readRDS("/projectnb2/talbot-lab-data/metabolic_models/scripts/N_cycle_sim/data/simulation_summary/more_sugars_total_biomass.rds")
biomass_long <- merge(biomass_long, eval_df[,c("sampleID","id")])







final_biomass = sim_biomass_wide %>% filter(cycle==max(cycle))

to_eval_abundances = merge(final_biomass, nitr_fg_abun, by = c("siteID","plot_date", "dateID", "plotID","horizon"))
to_eval_abundances = merge(to_eval_abundances, nitr_gen_abun, by = c("siteID","plot_date", "dateID", "plotID","horizon"))

to_eval_abundances %>%
	filter(cycle==300) %>%
	ggplot() + geom_point(aes(x = nitrifier_abundance, rel_NOB))

to_eval_abundances %>%
	ggplot() + geom_point(aes(x = nitrifier_abundance, rel_NOB, color = cycle))

to_eval_abundances %>%
	ggplot() + geom_point(aes(x = nitrifier_abundance, rel_AOB))

to_eval_abundances %>%
	ggplot() + geom_point(aes(x = nitrifier_abundance, rel_AOB_NOB))

to_eval_abundances %>%
	ggplot() + geom_point(aes(x = nitrifier_abundance, y= (Nitrospira_defluvii_NOB/all))

to_eval_abundances %>%
													ggplot() + geom_point(aes(x = nitrifier_abundance, y= rel_Nitrospira))



to_eval_metagenome = merge(final_biomass, nitrospiracae, by = c("siteID","plot_date", "dateID", "plotID","horizon"))

												to_eval_metagenome %>%
													ggplot() + geom_point(aes(x = Nitrospiraceae_pct, y= rel_Nitrospira))

												to_eval_metagenome %>%
													ggplot() + geom_point(aes(x = Nitrospiraceae_pct, y= rel_AOB_NOB))

												to_eval_abundances %>% ggplot(aes(x = horizon, y = rel_Nitrospira)) + geom_boxplot() + geom_jitter(alpha=.5)
												to_eval_abundances %>% ggplot(aes(x = horizon, y = rel_AOB_NOB)) + geom_boxplot() + geom_jitter(alpha=.5)
												to_eval_abundances %>% ggplot(aes(x = horizon, y = nitrifier_abundance)) + geom_boxplot() + geom_jitter(alpha=.5)

												to_eval_abundances %>% ggplot() + geom_boxplot(aes(x = horizon, y = rel_AOB_NOB))
												to_eval_abundances %>% ggplot() + geom_boxplot(aes(x = horizon, y = nitrifier_abundance))




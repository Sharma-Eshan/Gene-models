## code to read CIF files from AF database and calculate mean prediction score, total residues and and counts of plddt > 50, 60, 70, 80 and 90 for each structure. 

# In linux terminal , run the following to unzip .tar folder at UP000059680_39947_ORYSJ_v6.tar and unzip all files inside the folder
# curl -O https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000059680_39947_ORYSJ_v6.tar
# tar -xf ./UP000059680_39947_ORYSJ_v6.tar -C ./test_rice_af/
# cd ./test_rice_af
# gzip -d *.gz


# in Rstudio 

install.packages("bio3d")
install.packages("dplyr")
library(bio3d)
library(dplyr)

base_dir = "./test_rice_af"
pattern_cif <- ".cif$"

cif_files <- list.files(path = base_dir, pattern = pattern_cif, recursive = TRUE, full.names = TRUE)

cif_list <- list()

# This step may take time depending on the number of files
cif_list <- lapply(cif_files, function(f) {
  id <- sub(pattern_cif, "", basename(f))
  cif_list[[id]] <- read.cif(f)
  cif_list[[id]] <-  as.data.frame(cif_list[[id]]$atom)
  })


names(cif_list) <- sub(pattern_cif, "", basename(cif_files))

# from the b column of each df of the list cif_list, calculate mean value, number of total residues , number of residues with b > 50, 60, 70, 80 and 90

cif_stats <- lapply(names(cif_list), function(id) {
  df <- cif_list[[id]]
  df$b <- as.numeric(df$b)
  total_residues <- length(unique(df$resno))
  df1 <- df %>% group_by(resno) %>% summarise(plddt = mean(b)) %>% 
    summarise(mean_plddt = mean(plddt, na.rm = TRUE),
              plddt_50 = sum(plddt > 50),
              plddt_60 = sum(plddt > 60),
              plddt_70 = sum(plddt > 70),
              plddt_80 = sum(plddt > 80),
              plddt_90 = sum(plddt > 90))
  data.frame(
    id = id,
    mean_plddt = df1$mean_plddt,
    total_residues = total_residues,
    plddt_50 = df1$plddt_50,
    plddt_60 = df1$plddt_60,
    plddt_70 = df1$plddt_70,
    plddt_80 = df1$plddt_80,
    plddt_90 = df1$plddt_90,
    stringsAsFactors = FALSE
  )
})

cif_df <- do.call(rbind, cif_stats)

write.csv(cif_df, file = "Alphafold_stats_from_tar.csv", row.names = FALSE)

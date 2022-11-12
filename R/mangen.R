library(qqman)
##generate Manhattan plots for any traits which have any SNPS that pass correction
library(genio)
library(patchwork)
library(ggplotify)

figpath_individual <- "/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/gwas_results_figures/"
figpath_panels <- "/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/gwas_results_figures_panels/"
all_gwases <- c()
numGwases = length(all_gwases) ##correct for all the GWASES that we did, not just the ones from each loader.


##make all plots for traits with significant hits and put them in a directory
make_all_gwas_plots_individual <- function(clumped = TRUE){
  if (clumped){basepath <- "/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/gwas_results_clumped/"
  }
  else{
    basepath <- "/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/gwas_results/"

  }
    for (long_filename in list.files(basepath)){
    if (endsWith(long_filename, ".clumped"    )){

      if (clumped){gwas <- read_clumped(paste0(basepath, long_filename))}
      else{
        gwas <- read.csv(paste0(basepath, long_filename, sep = "\t"))
      }

      if (length(gwas$P) != 0){
        if (nrow(gwas[gwas$P < (5*10**(-8)),]) != 0){
          jpeg(stringr::str_replace_all(paste0(figpath_individual, long_filename, ".jpg"), "%", ","), width =1080, height =1080, quality =100)
          if (clumped){
            manhattan(gwas, chr = "CHR",
                      bp = "BP",
                      snp = "SNP",
                      yaxs = "i",
                      annotatePval = 5*10**(-8),
                      main = paste0("Manhattan Plot of ", stringr::str_replace_all(stringr::str_replace_all(long_filename, ".clumped", ""), "batch0.", "")),
                      col = c("navy", "steelblue"), ylim = c(4, max(-log10(gwas$P)))) ##the significance threshol dfor clumping is 0.0001, so we won't have any SNPS below this
          }
          else{
            manhattan(gwas, chr = "CHR",
                      bp = "BP",
                      snp = "SNP",
                      yaxs = "i",
                      annotatePval = 5*10**(-8),
                      main = paste0("Manhattan Plot of ", stringr::str_replace_all(stringr::str_replace_all(long_filename, ".clumped", ""), "batch0.", "")),
                      col = c("navy", "steelblue"))

          }
          dev.off()

    }

      }
    }
  }
}

##only for clumped traits now
panel_plots_panelby_loader<- function(loader = "percentile"){
  basepath <- "/net/mraid08/export/jasmine/zach/height_gwas/all_gwas/gwas_results_clumped/"
  library("patchwork")
  if (loader == "percentile"){
    operative_gwases <- c('batch0.in_range_70_180.clumped',
                         'batch0.MAD.clumped',
                         'batch0.iqr.clumped',
                         'batch0.below_70.clumped',
                         'batch0.sd.clumped',
                         'batch0.body_arm_left_bmd.clumped',
                         'batch0.grade_eugly.clumped',
                         'batch0.above_140.clumped',
                         'batch0.body_head_bmc.clumped',
                         'batch0.spine_l1_average_width.clumped',
                         'batch0.SdHHMM.clumped',
                         'batch0.spine_l1_l2_average_width.clumped',
                         'batch0.SdW.clumped',
                         'batch0.body_arms_bmd.clumped',
                         'batch0.body_arm_right_bmd.clumped',
                         'batch0.below_54.clumped')
  }
  else if (loader == "bloodtests"){
    operative_gwases <- c('batch0.bt__protein_total.glm.linear',
      'batch0.bt__bilirubin_total.glm.linear',
      'batch0.bt__uric_acid.glm.linear',
      'batch0.bt__basophils_%.glm.linear',
      'batch0.bt__glucose.glm.linear',
      'batch0.bt__monocytes_abs.glm.linear',
      'batch0.bt__vitamin_d.glm.linear',
      'batch0.bt__hdl_cholesterol.glm.linear',
      'batch0.bt__folic_acid.glm.linear',
      'batch0.bt__mean_platelet_volume.glm.linear',
      'batch0.bt__alkaline_phosphatase.glm.linear',
      'batch0.bt__tsh.glm.linear',
      'batch0.bt__monocytes_%.glm.linear')
  }
  else if (loader == "cgm"){
    operative_gwases <- c('batch0.Median.glm.linear',
      'batch0.lbgi.glm.linear',
      'batch0.j_index.glm.linear',
      'batch0.AUC.glm.linear',
      'batch0.adrr.glm.linear',
      'batch0.1st_Qu..glm.linear',
      'batch0.Min..glm.linear',
      'batch0.Mean.glm.linear',
      'batch0.3rd_Qu..glm.linear',
      'batch0.grade_hypo.glm.linear',
      'batch0.ea1c.glm.linear',
      'batch0.hypo_index.glm.linear',
      'batch0.gmi.glm.linear')
  }
  else if (loader == "ultrasound"){
    operative_gwases <- c('batch0.q_box_mean_mm_qbox_diameter.glm.linear',
      'batch0.q_box_mm_2_qbox_diameter.glm.linear',
      'batch0.q_box_median_mm_qbox_diameter.glm.linear',
      'batch0.q_box_mm_1_qbox_diameter.glm.linear',
      'batch0.q_box_mm_3_qbox_diameter.glm.linear')

  }
  else if (loader  == "metabolomics"){
    operative_gwases <- c('batch0.Lipids_POS_951.7732_412.1897_348.8896',
                           'batch0.Lipids_POS_772.5478_312.7962_299.3659',
                           'batch0.Lipids_POS_878.6730_357.8089_331.5334',
                           'batch0.Lipids_POS_925.7570_412.2867_343.8190',
                           'batch0.Lipids_POS_934.7356_388.2823_344.5864',
                           'batch0.Lipids_NEG_736.6819_427.4490_307.8527',
                           'batch0.Lipids_POS_587.5473_371.7741_270.3309',
                           'batch0.Lipids_POS_692.6319_447.7308_310.2593',
                           'batch0.Lipids_NEG_310.8379_48.6147_133.8410',
                           'batch0.Lipids_POS_644.6326_403.8752_292.8943',
                           'batch0.Lipids_NEG_720.6505_403.8259_300.8752',
                           'batch0.Lipids_POS_955.8045_433.6478_354.5816',
                           'batch0.Lipids_POS_722.7121_419.3648_314.3103',
                           'batch0.Lipids_NEG_1089.8691_431.3047_412.3494',
                           'batch0.Lipids_NEG_698.6211_419.3794_292.0184',
                           'batch0.Lipids_POS_953.7890_422.5644_352.0154',
                           'batch0.Lipids_NEG_984.8941_499.7480_364.6446',
                           'batch0.Lipids_POS_975.7739_407.4706_353.2466',
                           'batch0.Lipids_POS_927.7719_422.3191_347.1757')
  }
  ps <- list()
  i <- 1
  for (short_filename in operative_gwases){
    if (loader == "metabolomics"){
      long_filename <- paste0(short_filename, ".glm.linear")
    }
    else{
      ##for these loaders, the name is complete as is is
      long_filename <- short_filename
    }
      gwas <- read_clumped(paste0(basepath, long_filename))
      if (length(gwas$P) != 0){
          p <- ggplot_manhattan(gwas, theTitle = stringr::str_replace(short_filename, pattern = "batch.0", replacement = ""))
          ps[[i]] <- p
          i <- i + 1
      }
  }
  ##do this here
  jpeg(stringr::str_replace_all(paste0(figpath_panels, loader, ".jpg"), "%", ""), width = 1920, height = 3000)
  cowplot::plot_grid(plotlist = ps)
}



library(parallel)

# 模板文件路径
template_path <- "D:/XXYDATAanalysis/IP_1/ClusterFile_PC/code/mlGAR_modelfit_v2A_SingleStageNoCor_pc_template.R"

# 输出目录
output_dir <- "D:/XXYDATAanalysis/IP_1/ClusterFile_PC/results"

# 参数设置
params <- list(
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE = 1),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE = 2),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE = 3),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE = 4),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =5),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =6),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =7),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =8),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =9),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =10),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =11),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =12),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =13),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =14),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =15),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =16),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =17),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =18),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =19),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =20),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =21),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =22),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =23),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =24),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =25),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =26),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =27),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =28),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =29),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =30),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =31),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =32),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =33),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =34),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =35),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =36),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =37),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =38),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =39),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =40),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =41),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =42),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =43),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =44),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =45),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =46),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =47),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =48),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =49),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =50),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =51),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =52),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =53),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =54),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =55),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =56),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =57),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =58),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =59),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =60),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =61),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =62),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =63),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =64),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =65),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =66),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =67),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =68),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =69),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =70),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =71),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =72),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =73),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =74),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =75),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =76),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =77),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =78),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =79),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =80),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =81),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =82),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =83),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =84),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =85),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =86),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =87),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =88),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =89),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =90),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =91),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =92),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =93),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =94),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =95),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =96),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =97),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =98),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =99),
  list(nT_REPLACE=15, nP_REPLACE = 150, r_REPLACE =100),
  
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE = 1),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE = 2),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE = 3),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE = 4),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =5),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =6),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =7),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =8),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =9),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =10),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =11),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =12),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =13),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =14),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =15),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =16),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =17),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =18),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =19),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =20),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =21),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =22),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =23),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =24),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =25),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =26),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =27),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =28),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =29),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =30),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =31),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =32),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =33),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =34),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =35),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =36),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =37),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =38),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =39),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =40),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =41),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =42),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =43),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =44),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =45),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =46),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =47),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =48),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =49),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =50),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =51),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =52),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =53),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =54),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =55),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =56),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =57),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =58),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =59),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =60),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =61),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =62),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =63),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =64),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =65),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =66),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =67),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =68),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =69),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =70),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =71),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =72),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =73),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =74),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =75),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =76),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =77),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =78),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =79),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =80),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =81),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =82),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =83),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =84),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =85),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =86),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =87),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =88),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =89),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =90),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =91),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =92),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =93),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =94),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =95),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =96),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =97),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =98),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =99),
  list(nT_REPLACE=50, nP_REPLACE = 500, r_REPLACE =100)
  
)

# 读取模板文件
template <- readLines(template_path)

# 函数：创建并运行脚本
create_and_run_script <- function(params_temp) {#id=1
  new_script <- template
  
  # 替换变量
  for (param_name in names(params_temp)) {
    new_script <- gsub(paste0("\\b", param_name, "\\b"), params_temp[[param_name]], new_script)
  }
  
  # 构建文件名
  nP <- params_temp$nP_REPLACE
  nT <- params_temp$nT_REPLACE
  r <- params_temp$r_REPLACE
  # 定制文件名
  script_name <- paste0("mlGAR_v2A_SingleStageNoCor_nP", nP, "nT", nT,"r",r,".R")
  script_path <- file.path(output_dir, script_name)
  
  # 保存新脚本
  writeLines(new_script, script_path)
  
  # 运行新脚本
  message("Running: ", script_path)
  system(paste("Rscript", script_path))
}

# 设置并行核心数
numCores <- 7 #detectCores() - 1

# 使用parLapply并行运行所有脚本
cl <- makeCluster(numCores)
clusterExport(cl, list("create_and_run_script", "template", "output_dir", "params"))

results <- parLapply(cl, params, create_and_run_script)
stopCluster(cl)

# 查看结果
print(results)
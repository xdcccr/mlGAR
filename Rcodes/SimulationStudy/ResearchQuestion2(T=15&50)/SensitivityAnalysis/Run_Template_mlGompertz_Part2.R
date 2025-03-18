library(parallel)

# 模板文件路径
template_path <- "D:/XXYDATAanalysis/IP_1/resultsRevision1/1MLGompertzWithMoreInformativePriors/Template_mlGompertz_Part2.R"

# 输出目录
output_dir <- "D:/XXYDATAanalysis/IP_1/resultsRevision1/results_1"

# 参数设置
params <- list(
  list(nT_REPLACE=15, nP_REPLACE = 150),
  
  list(nT_REPLACE=50, nP_REPLACE = 500)
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
  # Round <- params_temp$Round_REPLACE
  
  # 定制文件名
  # script_name <- paste0("mlGompertz_nP", nP, "nT", nT,"Round",Round,".R")
  script_name <- paste0("mlGompertz_nP", nP, "nT", nT,".R")
  script_path <- file.path(output_dir, script_name)
  
  # 保存新脚本
  writeLines(new_script, script_path)
  
  # 运行新脚本
  message("Running: ", script_path)
  system(paste("Rscript", script_path))
}

# 设置并行核心数
numCores <- detectCores() - 1

# 使用parLapply并行运行所有脚本
cl <- makeCluster(numCores)
clusterExport(cl, list("create_and_run_script", "template", "output_dir", "params"))

results <- parLapply(cl, params, create_and_run_script)
stopCluster(cl)

# 查看结果
print(results)
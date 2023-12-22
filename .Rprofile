vscode_init_path <- file.path(Sys.getenv(
    if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"
), ".vscode-R", "init.R")
if (file.exists(vscode_init_path)) {
    source(vscode_init_path)
}
# /Users/joankant/Desktop/gaitigroup/Users/Joan/h4h-cell-cell-interactions/.Rprofile
if (Sys.info()[["sysname"]] == "Darwin") {
    # If for mac store on local machine
    Sys.setenv(RENV_PATHS_CACHE = "~/Library/CloudStorage/OneDrive-UHN/renv_cache")
    Sys.setenv(CLANG = "~/miniforge3/envs/standard_env/bin/clang")
    Sys.setenv("CLANG++" = "~/miniforge3/bin/clang++")
    Sys.setenv("pkg-config" = "~/miniforge3/envs/standard_env/bin/pkg-config")
    Sys.setenv(GCC = "~/miniforge3/envs/standard_env/bin/gcc")
    Sys.setenv("libxext" = "~/miniforge3/envs/standard_env/lib/libXext.dylib")
    Sys.setenv(RENV_PATHS_LIBRARY = "renv/library")
} else if (Sys.info()[["sysname"]] == "Linux") {
    # If for linux store on cluster
    # if(dir.exists)
    # TODO: change to project directory
    project_dir <- "h4h-cell-cell-interactions"
    if (!dir.exists(paste0("~/renv_libs/", project_dir, "_library"))) {
        dir.create(paste0("~/renv_libs/", project_dir, "_library"), recursive = TRUE)
    }
    Sys.setenv(RENV_PATHS_CACHE = paste0("~/renv_cache"))
    Sys.setenv(RENV_PATHS_LIBRARY = paste0("~/renv_libs/", project_dir, "_library"))
}

# TODO: testing
.libPaths(c(
    .libPaths(), "~/Desktop/gaitigroup/Users/Joan/h4h-cell-cell-interactions/renv/library",
    "~/miniforge3/envs/standard_env/lib/R/library"
))

options(
    renv.config.sandbox.enabled = FALSE,
    renv.config.cache.enabled = TRUE,
    renv.settings.use.cache = TRUE,
    renv.consent = TRUE,
    renv.config.auto.snapshot = FALSE,
    renv.config.pak.enabled = FALSE
)
# source("renv/activate.R")

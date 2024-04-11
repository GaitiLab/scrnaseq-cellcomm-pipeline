vscode_init_path <- file.path(Sys.getenv(
    if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"
), ".vscode-R", "init.R")
if (file.exists(vscode_init_path)) {
    source(vscode_init_path)
}

options(future.globals.maxSize = 8000 * 1024**2)

if (interactive()) {
    suppressMessages(require(devtools))
}

suppressMessages(library(triplediff))
df <- read.csv("C:/Users/35223/Desktop/默认项目/程序包/ddd/ddd_sim.csv")
df$id <- as.numeric(df$id)

run <- function(cg, xf) {
  res <- tryCatch(
    triplediff::ddd(yname="y", tname="time", idname="id", gname="enable",
                    pname="eligible", xformla=xf, data=df,
                    control_group=cg, panel=TRUE),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(res)) {
    cat("control_group =", cg, " xformla =", deparse(xf), "\n")
    print(summary(res))
  }
}

run("nevertreated", ~1)

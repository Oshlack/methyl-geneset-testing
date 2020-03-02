library(workflowr)
library(here)
wflow_publish(files = here("analysis/regionAnalysis.Rmd"))

paletteer_display <- function(pkg){

       colorlist <- palettes_d_names$palette[palettes_d_names$package == pkg]
       nr <- length(colorlist)
       nc <- max(palettes_d_names$length[palettes_d_names$package == pkg])
       n <- palettes_d_names$length[palettes_d_names$package == pkg]
       ylim <- c(0, nr)
       oldpar <- par(mgp = c(2, 0.25, 0))
       on.exit(par(oldpar))
       plot(1, 1, xlim = c(0, nc), ylim = ylim, type = "n", axes = FALSE,
                    bty = "n", xlab = "", ylab = "")
       for (i in 1:nr) {
           nj <- n[i]
           if (colorlist[i] == "")
               next
           shadi <- paletteer_d(package = !!pkg, palette = !!colorlist[i])
           rect(xleft = 0:(nj - 1), ybottom = i - 1, xright = 1:nj,
                          ytop = i - 0.2, col = shadi, border = "light grey")
         }
       text(rep(-0.1, nr), (1:nr) - 0.6, labels = colorlist, xpd = TRUE,
                    adj = 1)
     }

pdf(here("output/paletteer-palettes.pdf"))
sapply(unique(palettes_d_names$package), function(p){
  paletteer_display(p)
})
dev.off()

setwd("~/Documents/IC/Lipidome/Data")

library(LOBSTAHS)
library(CAMERA)
library(xcms)

xs <- xcmsSet(files = "BD01.mzXML", method = "centWave")
xsa <- xsAnnotate(xs)
anF <- groupFWHM(xsa)
anI <- findIsotopes(anF, ppm = 2)
rm(xs, xsa, anF)
anIC <- groupCorr(anI)
anFA <- findAdducts(anIC, polarity = "negative")

data <- doLOBscreen(anFA,
                    polarity = "negative",
                    match.ppm = 2)

peaklist <- peakdata(data)
write.csv(peaklist, file = "xsannotated.csv")


library(IPO)
peakpickingParameters = getDefaultXcmsSetStartingParams('centWave')
resultPeakpicking = optimizeXcmsSet(files = "BD01.mzXML",
                                    params = peakpickingParameters,
                                    nSlaves=4,
                                    subdir='rsmDirectory')

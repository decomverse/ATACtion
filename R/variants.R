load_gwascat_GR <- function(current = TRUE, coordinate = 'hg38') {	
	library(gwascat)
	
	if(current == TRUE)
		suppressWarnings({ebicat38 = gwascat::makeCurrentGwascat()})
	else
		data("ebicat38", package = "gwascat", envir = environment())
	
	if(coordinate == 'hg19') { # Liftover
		ebicat19 = liftOverGR(ebicat38, "hg38", "hg19")
		return(as(ebicat19, "GRanges"))
	} else {
		return(as(ebicat38, "GRanges"))
	}
}

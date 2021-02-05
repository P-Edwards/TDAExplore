Rcpp::loadModule("Landscape", TRUE)

#' Compute landscape of a persistence diagram.
#'
#' @param PersistenceDiagram A pd object from the diagram function, or a list of pairs.
#' @param exact Set to TRUE for exact computation, FALSE for discrete.
#' @param max_x Domain threshold for discrete PL.
#' @param dx Domain grid diameter for discrete PL.
#' @param degree If input is pd object then degree specifies which degree to select from.
#' @param threshold Threshold used to compute PD (could be inf).
#' @return A PersistenceLandscape object.
landscape <- function(PersistenceDiagram, degree=NULL, exact=FALSE, dx=0.1,  min_x=0, max_x=10, threshold=-1){

    diagram = NULL
    max_y = threshold
    
    #Automatic parameter deduction logic:
    
    # Input is bd pairs (no threshold info).
    if( is.atomic(PersistenceDiagram) || is.null(PersistenceDiagram$param.threshold)){
        diagram = PersistenceDiagram

        if(threshold == -1){
            max_y = max(diagram[,2])
        }
    }
    
    # Input is from diagram output.
    else{
        if(is.null(degree)){
            stop("Error: If input persistence diagram is directly from 
                 the diagram fucntion then a homological degree degree must be specified (add degree=__ to landscape call). ")
            
        }

        else{
            diagram = PersistenceDiagram$pairs[[degree]]
            max_y = PersistenceDiagram$param.threshold
        }
    }

    
    #sanity checks
    if(is.null(diagram) || all(is.na(diagram))){
        stop("Error: Empty persistence diagram.")
    }

	#Construct a persistence landscape.
    landscape_raw <- methods::new(PersistenceLandscape, diagram, exact, min_x, max_x, dx, max_y)

    return(landscape_raw)
}



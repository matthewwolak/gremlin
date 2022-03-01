#' Optimization Algorithm Checks.
#'
#' Check and/or set optimization algorithms to use. Intended for internal use
#' within gremlin
#'
#' @aliases algChk
#' @param algit A character \code{vector} of algorithms for each iteration.
#' @param maxit An \code{integer} indicating the number of REML iterations.
#' @param ctrl A \code{list} of arguments set by \code{gremlinControl}.
#' @param mc A model \code{call}.
#'
#' @return A \code{list} containing a vector specifying the currently implemented
#'   optimization algorithms for each iteration and a vector containing the type
#'   of finite differences for each iteration if this procedure is to be used for
#'   first derivative calculations.
#' @author \email{matthewwolak@@gmail.com}
algChk  <- function(algit, maxit, ctrl, mc){
  # List of all acceptable algorithm choices
  #algChoices <- c("EM", "AI", "AIbfd", "AIcfd", "AIffd", "bobyqa", "NR", ctrl$algorithm)
  algChoices <- c("EM", "AI", "AIbfd", "AIcfd", "AIffd", "bobyqa", "NR")  #<-- ignore ctrl$algorithm for now

    if(!is.null(ctrl$algorithm)){
      #TODO check validity of `ctrl$algorithm` and `ctrl$algArgs`
      ## need to pass algorithm to `gremlinR` or switch to it if `gremlin` called
      ## temporarily IGNORE with warning
      warning(cat("Ignored algorithm(s) supplied in control",
        dQuote(ctrl$algorithm),
        ". gremlin is not old enough for user-specified algorithms\n"),
          immediate. = TRUE)
    }
    
  algMatch <- pmatch(algit, algChoices, nomatch = 0, duplicates.ok = TRUE)
  if(all(algMatch == 0) & !all(algit %in% ctrl$algorithm)){ 
      stop(cat("Algorithms:", dQuote(algit[which(algMatch == 0)]),
        "not valid. Please check values given to the `algit` argument\n"))
  }
  if(any(algMatch == 0)){  
    warning(cat("Algorithms:", dQuote(algit[which(algMatch == 0)]),
      "not valid - dropped from the list\n"))
    algit <- algit[-which(algMatch == 0)]
  }
  if(is.null(mc$algit)){
#    algit <- c(rep("EM", min(maxit, 2)), rep("AIcfd", max(0, maxit-2)))
    algit <- c("AIcfd")
  } else algit <- algChoices[algMatch]
#  if(length(algit) == 0) algit <- c(rep("EM", min(maxit, 2)),
#                                    rep("AIcfd", max(0, maxit-2)))
  if(length(algit) == 0) algit <- c("AIcfd")
  if(length(algit) == 1) algit <- rep(algit, maxit)
  # Now check for finite difference algorithms with AI:
  fdit <- as.integer(rep(0, maxit))
    fdit[grep("bfd", algit)] <- 1
      algit <- gsub("bfd", "fd", algit)
    fdit[grep("cfd", algit)] <- 2
      algit <- gsub("cfd", "fd", algit)
    fdit[grep("ffd", algit)] <- 3
      algit <- gsub("ffd", "fd", algit)
      
 return(list(algit = algit, fdit = fdit))
} 


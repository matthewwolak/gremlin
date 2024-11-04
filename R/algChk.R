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
#' @return A \code{list} containing a vector \code{algit} specifying the currently
#'   implemented optimization algorithms for each iteration along with a vector
#'   each containing the type of first (\code{fdit}) and second (\code{sdit})
#'   derivatives for each iteration (else \code{NA} if either is not applicable).
#' @author \email{matthewwolak@@gmail.com}
#' @export
algChk  <- function(algit, maxit, ctrl, mc){
  # List of all acceptable algorithm choices
  #algChoices <- c("EM", "AI", "AItr", "AIbfd", "AIcfd", "AIffd",
  #  "bobyqa", "NR", ctrl$algorithm)
  algChoices <- c("EM", "AI", "AItr", "AIbfd", "AIcfd", "AIffd", "bobyqa", "NR")  #<-- ignore ctrl$algorithm as last entry for now

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
        "not valid. Check values given to the `algit` argument are one of:",
        algChoices, "\n"))
  }
  if(any(algMatch == 0)){  
    warning(cat("Algorithms:", dQuote(algit[which(algMatch == 0)]),
      "not valid - dropped from the list\n"))
    algit <- algit[-which(algMatch == 0)]
  }
  if(is.null(mc$algit) & is.null(algit)){
    algit <- c(rep("EM", min(maxit, 1)), rep("AItr", max(0, maxit-1)))
  } else algit <- algChoices[algMatch]
  if(length(algit) == 0) algit <- c(rep("EM", min(maxit, 1)),
                                    rep("AItr", max(0, maxit-1)))
  # Swap out "AI" for explicit/default first derivative method
  if(length(AI <- grep("(^|\\s)AI($|\\s)", algit, ignore.case = FALSE)) > 0){
    algit[AI] <- "AItr"  
  }
  if(length(algit) < maxit) algit <- rep_len(algit, length.out = maxit)
  
  # Now check for first and second derivative calculations/algorithms:
  fdit <- sdit <- rep(NA, maxit)
    fdit[grep("bfd", algit)] <- "bfd"   #1
    fdit[grep("cfd", algit)] <- "cfd"   #2
    fdit[grep("ffd", algit)] <- "ffd"   #3
    fdit[grep("tr", algit)] <- "tr"     #4
    fdit[grep("EM", algit)] <- "dfree"  #5

  sdit[grep("AI", algit)] <- "AI"
  #TODO something for either ctrl$alg or bobyqa/NR/etc.
          
 return(list(algit = algit,
   fdit = factor(fdit,
     levels = c("bfd", "cfd", "ffd", "tr", "dfree"), ordered = TRUE),
   sdit = sdit))
} 


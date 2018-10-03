#FIXME: Make decision about how to convert into factors (do this in function passing data frame to `mkModMats()`).
#FIXME: previous to `mkModMats()` NEED to make all rand effects/intercepts factors in data!
## Somehow do this in the data.frame associated with the model.frame in gf?

## For example, leave variables 'as is' in fixed effects part
## Convert all variables in certain parts of random effects formulas to factors
### Apparently, lme4 does this automatically - see R-sig-ME thread 5 October 2016 
#### Subject: manage person identifier variable
#### XXX Note the caution about leading zeroes! (can `as.factor(as.character())` solve this?)
## For example, year in both fxd and rfx (covariate as fixed effect, but factor as rfx)


####	start function		########
  # Strip out function so most of checking and cleaning done elsewhere
  # see "./gremlin_parent/modelFittingInR" for general strategy of model.frame and model.matrix
  # check for:
    # multivariate/univariate
    # sorting of data frame (multivariate or by interactions of random effects)   
    # if levels of fixed/random effects in their interaction are all present (e.g., MCMCglmm does this)
  #Have "units" added to data.frame in a previous function so it is already in "data="object
  #Perhaps remove some of arguments below that are more checking, constructing data

#TODO: `mkModMats()` that accepts lme4 model terms (e.g., "|" for rand fx)
# Below is the MCMCglmm kind of way




#XXX RANDOM EFFECTS MUST ALREADY BE FACTORS IN DATA!!!!
################################################################################
#' @describeIn gremlinR Generates model matrices. 
#' @export
#' @importFrom stats na.pass model.response model.matrix lm reformulate na.omit
mkModMats <- function(formula, random = NULL, rcov = ~ units,
		data = NULL, subset = NULL,
		ginverse = NULL,
		na.action = na.pass, offset = NULL, contrasts = NULL,
		Xsparse = TRUE,   # Should fxd effects design matrix be sparse
		...)
{

  cl <- match.call()
  ## keep only the arguments which should go into the model frame
  #TODO: Need to build up 1 model frame with all of the components for the fit
  ## perhaps can add back X, y, Z, etc to 'mf' and return 'mf' - as list of fixed and random frames?

  if(any(grepl("units", names(data)))){
    stop("column names in data cannot be named 'units'")
  }
  if(!any("units" %in% all.vars(rcov))){
    stop("'units' must be specified in 'rcov' argument")
  }
  #FIXME: remove next stop when R structures can be specified
  if(deparse(rcov) != "~units"){
    stop("'rcov' can only specify a single R variance at this time")
  }




  ##################################
  #####		Fixed effects  	####
  mf <- match.call(expand.dots = FALSE)
  f <- match(c("formula", "data", "subset", "na.action", "offset"), names(mf), 0)
  ff <- mf[c(1, f)]
  #TODO: How does this handle missing data in fixed effects?
  ff[[1]] <- as.name("model.frame")
  rf <- gf <- ff    			# create rand effects model frames 
  ff$drop.unused.levels <- TRUE         # do after assign gf so this is not set for rfx
  ff <- eval.parent(ff)   	        # fxd effects model frame

  ## 1) allow model.frame to update the terms object before saving it.
  ft <- attr(ff, "terms") 
  resp <- model.response(ff, "numeric")
    if(is.null(ncol(resp))) ncy <- 1 else ncy <- ncol(resp)  #<--number of traits
if(ncy > 1) stop("gremlin isn't old enough to play with multivariate models") #FIXME!!!
  #TODO: if multivar; order y as traits within individs (Meyer 1991, just b4 eqn 5)
  y <- matrix(resp, ncol = 1)
  if(ncy == 1) ny <- nrow(y) else ny <- NULL #FIXME  # <-- sample size if univariate

  ## 2) retrieve the weights and offset from the model frame so
  ## they can be functions of columns in arg data.
  #XXX offset <- model.offset(ff)                           # turned off for now

  #TODO: Need to make list (e.g. rand fx) and create X for each response (in multivar)
  ## bdiag(lapply(seq(ncol(resp)), FUN = function(x){model.matrix(...
  #TODO: Test if sparse or dense X is better/more efficient
  if(Xsparse){
    X <- sparse.model.matrix(ft, ff,
	contrasts,
	transpose = FALSE)
  } else X <- model.matrix(ft, ff,
	contrasts)
  ## if any subsetting is done, retrieve the "contrasts" attribute here.
  #TODO: Figure out what this means in original context of lm and if applies
  #TODO: use lm to find singularities in X (see `MCMCglmm()`, ~line 608 of v2.21)
### FIXME ##
sing.rm<-lm(y~ as(X, "matrix")-1)
sing.rm<-which(is.na(sing.rm$coef))
if(length(sing.rm)){
   warning("Some fixed effects are not estimable and have been removed/set to 0")
   X<-X[,-sing.rm, drop = FALSE]
}
###END FIXME
  #TODO OR created reduced rank X (see "./Wellham2008REML..." slide 15&16)
  #TODO see how lme4 deals with rank-deficiencies:
  ## `rankMatrix()` test function and discussion https://github.com/lme4/lme4/pull/144





  #####################################
  ######	R structure	#######
  
  #TODO: strip R formula of covariance functions (e.g., idh(), us(), etc.)
  rForm <- rcov
  #TODO: check below that when complex R structure this gets the formula right
  rf$formula <- reformulate(deparse(rForm[[-1]]), cl$formula[[2]], intercept = FALSE)
  #TODO: Need to add "units" to 'data'
  Rdata <- data
  Rdata$units <- factor(paste("units", seq(nrow(Rdata)), sep = "."), levels = paste("units", seq(nrow(Rdata)), sep = ".")) 
  rf$data <- Rdata
  rf$drop.unused.levels <- TRUE
  #TODO: how handle missing values in random effects?
  rf <- eval.parent(rf)   	# R model frame
  rt <- attr(rf, "terms") 
  Zr <- sparse.model.matrix(rt, rf,
	transpose = FALSE,
	row.names = TRUE) 



  ######################################
  ####		G structure	  ######
  if(!is.null(random)){
    if(random[[1]] != "~") stop("Tilde not first element in random formula")
  #TODO: Make a more rigorous check of formula (and if a formula == is.formula)



  # Split random formula up by terms "+" and pass, one by one to varStructures
#XXX: See what `terms(random)` does, particularly in a formula with multiple `us(trait):`
#XXX See how lme4 looks for terms in the formula, see `nobars`, but in full context of code
## what about `all.vars()` or `all.names()` for above questions?
#XXX (commented out 18 Oct 2015) varStructures(random[[2]])



  # Check the ginverse
  ## taken from MCMCglmm.R
    if(!is.null(ginverse)){
      if(!is.list(ginverse)){
        stop("ginverse must be an object of class 'list'")
      }
      for(i in 1:length(ginverse)){           
        if(is.null(rownames(ginverse[[i]]))){
          stop(paste(names(ginverse)[i], "ginverse must have non-null rownames"))
        }
        if(names(ginverse)[i] %in% names(data) == FALSE){
          stop(paste(names(ginverse)[i], "does not appear in data"))
        }
        dat_ilevels <- na.omit(unique(data[,names(ginverse)[i]]))
        if(any(dat_ilevels %in% rownames(ginverse[[i]]) == FALSE)){
          stop(cat("levels of", names(ginverse)[i], "do not have a row entry in ginverse:", as.character(dat_ilevels[dat_ilevels %in% rownames(ginverse[[i]]) == FALSE]), "\n"))
        }
        if(any(duplicated(rownames(ginverse[[i]])))){
          stop(paste("rownames of ", names(ginverse)[i], "ginverse must be unique"))
        }
#        if(determinant(ginverse[[i]])$sign<=0){
#          stop(paste(names(ginverse)[i], "ginverse is not positive definite"))
#        }
      #TODO: Now, check/turn ginverse list into sparse matrix (so that can allow ginverse elements to be triplet form
        #XXX For now, make sure ginverse is "dgCMatrix" - consider dsCMatrix to save space?
        if(class(ginverse[[i]]) != "dgCMatrix"){
          stop(paste(names(ginverse)[i], "ginverse of", class(ginverse[[i]]), "must be a 'dgCMatrix' class"))
        }
      ## Need to decide on "dgCMatrix" or "dsCMatrix"...
        #FIXME: Make correspond to how previous function makes factors/checks levels
        ## Need to change data either in call/function environment of the random effects model frame (gf)
        #below line: TEMPORARY to get factor levels in correct order so Z formed correctly
        data[,names(ginverse)[i]] <- factor(data[,names(ginverse)[i]], levels=rownames(ginverse[[i]]))                      
      }
    }
    #####


  # Make list of G and R elements as well as list of Zs (just for G effects?).
  ## Construction of MMEs or likelihood function can go element by element in a list as the rows and columns


## May need/be able to get rid of much of what is below and the model frame parts...
    #TODO: strip G formula of covariance functions (e.g., idh(), us(), etc.)
    #TODO: how to deal with random regression
    gForm <- random
    gf$formula <- gForm
    #TODO: how handle missing values in random effects?
    # do I need the response variable so that the random effects below get cut to only
    ## the effects that have a response variable??? 21 Oct 2015, following Mrode ch. 3 YES
    ### INCLUDE levels in Z that do not have y observation as columns, but not rows
    ### see how asreml's na.omit.Y creates Z...
    ###XXX See MCMCglmm's `buildZ()`, line 115 of v2.21 - drop levels unless Va term.
    #### Should I extend to any ginverse term (so V_dominance, etc.)
    #TODO: See MCMCglmm's `buildZ()` - why create t(Z)Z directly??
    Gnames <- strsplit(deparse(gForm[[-1]]), split = " \\+ ")[[1]]
    nG <- length(Gnames)
    #TODO: logDetG as a simple vector or scalar (1 element) that is added to!
    Zg <- listGeninv <- logDetG <- vector("list", length = nG)
    names(Zg) <- names(listGeninv) <- names(logDetG) <- Gnames

    for(g in 1:nG){
      ggf <- gf # temporary gf
      #XXX Does including response in `reformulate()` work when multivariate?
      ggf$formula <- reformulate(names(Zg)[g], response = cl$formula[[2]], intercept = FALSE)
      # retain levels in ginverse if one is associated, else drop unused levels
      if(!is.null(ginverse) && Gnames[g] %in% names(ginverse)){
        Ggdata <- data
        levels(Ggdata[, Gnames[g]]) <- union(ginverse[[g]]@Dimnames[[1]], levels(Ggdata[, Gnames[g]])) 
        ggf$data <- Ggdata
        ggf$drop.unused.levels <- FALSE
        ggf <- eval.parent(ggf)	# G[g] model frame
        gt <- attr(ggf, "terms")
        Zg[[g]] <- sparse.model.matrix(gt, ggf,
		transpose = FALSE,
		drop.unused.levels = FALSE,	#TODO correct?
		row.names = TRUE)
        listGeninv[[g]] <- ginverse[[which(names(ginverse) == Gnames[g])]]
        # calculate log(det(G)) from geninv; log(det(G)) = -1*log(det(G^-1))
	#TODO: make check to see if ginverse[[g]] has attribute=logdet (from nadiv)
        logDetG[[g]] <- -1 * determinant(ginverse[[which(names(ginverse) == Gnames[g])]], logarithm = TRUE)$modulus[1]
      } else{
          Zg[[g]] <- sparse.model.matrix(ggf$formula, eval.parent(ggf),
		transpose = FALSE,
		drop.unused.levels = TRUE,	#TODO correct?
		row.names = TRUE)
          #XXX Might want to do this in c++, see cs_kroneckerI
          #TODO check `n = ` value!?
          ## listGeninv need to be same dimensions as crossprod (in W) when making MMA
          listGeninv[[g]] <- Diagonal(x = 1, n = ncol(Zg[[g]])) 
          logDetG[[g]] <- 0
        }

    }


#FIXME better check for missing values or handling of missing values
## e.g., What does MCMCglmm do?
nrowZi <- sapply(seq(length(Zg)), FUN = function(g){nrow(Zg[[g]])})
if(any(nrowZi != ny)){
  missingZi <- which(nrowZi != ny)
  stop(cat("Missing value(s) for", names(Zg)[missingZi], "\n"))
}


  } else{
      Zg <- listGeninv <- logDetG <- NULL
      nG <- 0
    }


 structure(list(y = y, ny = ny, ncy = ncy,
	X = X, nb = ncol(X),
	Zr = Zr,
	Zg = Zg, nG = nG, listGeninv = listGeninv, logDetG = logDetG),
	class = "gremlinModMats")
}



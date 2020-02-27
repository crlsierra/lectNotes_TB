modWengLuo<-function
  ### This function creates an ecosystem model as described by Weng & Luo 2011.
  (t,      ##<< A vector containing the points in time where the solution is sought.
   U,	##<< Photosynthetically fixed carbon.
   b, ##<< Vector of partitioning coefficients of the photosynthetically fixed carbon.
   A, ##<< Matrix of transfer coefficients.
   C,	##<< Diagonal matrix with transfer coefficients.
   X0,    ##<< Vector with initial conditions.
   xi=1,  ##<< A scalar or data.frame object specifying the external (environmental and/or edaphic) effects on decomposition rates.
   solver=deSolve.lsoda.wrapper,  ##<< A function that solves the system of ODEs. This can be \code{\link{euler}} or \code{\link{ode}} or any other user provided function with the same interface.
   pass=FALSE  ##<< if TRUE Forces the constructor to create the model even if it is invalid 
  )	
  { 
    t_start=min(t)
    t_end=max(t)
    if(length(A)!=64) stop("A must be a matrix of dimention 8 x 8")
    if(length(C)!=64) stop("C must be a matrix of dimention 8 x 8")
    if(length(X0)!=8) stop("the vector with initial conditions must be of length = 8")
    if(length(b)!=8) stop("b must be of length = 8")
    
    if(length(U)==1){
      inputFluxes=BoundInFluxes(
        function(t){matrix(nrow=8,ncol=1,U*b)},
        t_start,
        t_end        
      )
    }
    if(class(U)=="data.frame"){
      x=U[,1]  
      y=U[,2]  
      inputFlux=splinefun(x,y)
      inputFluxes=BoundInFlux(
        function(t){matrix(nrow=8,ncol=1,inputFlux(t)*b)},
        min(x),
        max(x)        
      )
    }
    
    
    if(length(xi)==1) fX=function(t){xi}
    if(class(xi)=="data.frame"){
      X=xi[,1]
      Y=xi[,2]
      fX=splinefun(X,Y)
    }
    Af=BoundLinDecompOp(
      function(t){fX(t)*(A%*%C)},
      t_start,
      t_end      
    )
    
    Mod=GeneralModel(t=t,A=Af,ivList=X0,inputFluxes=inputFluxes,pass=pass)
    return(Mod)
    ### A Model Object that can be further queried 

  }
  

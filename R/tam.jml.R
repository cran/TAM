

tam.jml <- function( resp , group = NULL , adj=.3 , disattenuate = FALSE ,
            bias = TRUE, xsi.fixed=NULL ,  xsi.inits = NULL ,  theta.fixed = NULL , 
            A=NULL , B=NULL , Q=NULL , ndim=1 ,
            pweights = NULL , control = list() , version = 2 )
{           
	#**** version = 1
	if (version == 1){
		res <- tam_jml_version1( resp=resp , group = group , adj=adj , 
				disattenuate = disattenuate ,
				bias = bias, xsi.fixed=xsi.fixed ,  xsi.inits = xsi.inits ,  
				A=A , B=B , Q=Q , ndim=ndim , theta.fixed = theta.fixed , 
				pweights = pweights , control = control  )
	}
	#**** version = 2
	if (version == 2){
		res <- tam_jml_version2( resp=resp , group = group , adj=adj , 
				disattenuate = disattenuate ,
				bias = bias, xsi.fixed=xsi.fixed ,  xsi.inits = xsi.inits ,  
				A=A , B=B , Q=Q , ndim=ndim ,
				pweights = pweights , control = control  )
	}
	return(res)	
}
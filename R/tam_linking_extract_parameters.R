
tam_linking_extract_parameters <- function( tamobj )
{
	A <- tamobj$A
	xsi <- tamobj$xsi$xsi
	guess <- tamobj$guess	
	AXsi <- tamobj$AXsi
	B <- tamobj$B
	dimnames(A)[[1]] <- rownames(AXsi) <- dimnames(B)[[1]]		
	ndim <- dim(B)[3]
	items <- colnames(tamobj$resp)
	M <- tamobj$beta[1,1]
	SD <- sqrt( tamobj$variance[1,1] )	
	class_tamobj <- class(tamobj)
	#--- OUTPUT
	res <- list(A=A, xsi=xsi, guess=guess, AXsi=AXsi, B=B, ndim=ndim, items=items, M=M, SD=SD ,
					class_tamobj=class_tamobj)
	return(res)
}
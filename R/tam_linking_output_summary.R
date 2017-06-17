
tam_linking_output_summary <- function( parameters_list, linking_list , NM )
{
	#--- means and SDs
	M_SD <- matrix( 0 , nrow=NM , ncol=2)
	colnames(M_SD) <- c("M", "SD")
	rownames(M_SD) <- paste0( "study",1:NM)
	M_SD <- as.data.frame(M_SD)
	for (mm in 1:NM){
		M_SD[mm, "M"] <- parameters_list[[mm]]$M
		M_SD[mm, "SD"] <- parameters_list[[mm]]$SD		
	}
	#--- transformation item parameters
	trafo_items <- 0*M_SD
	colnames(trafo_items) <- c("a","b")
	trafo_items$a <- 1
	for (mm in 1:(NM-1) ){
		trafo_items[mm+1,] <- linking_list[[mm]]$linking_results$trafo_items
	}
	#--- transformation person parameters
	trafo_persons <- trafo_items
	for (mm in 1:(NM-1) ){
		trafo_persons[mm+1,] <- linking_list[[mm]]$linking_results$trafo_persons
	}	
	#--- number of linking items
	N_common <- matrix(0 , nrow=NM-1, ncol=1)
	rownames(N_common) <- paste0("Linking Study " , 1:(NM-1) , " -> Study " , 2:NM )
	colnames(N_common) <- "N_Items"
	for (mm in 1:(NM-1) ){
		N_common[mm,1] <- length( linking_list[[mm]]$common_items )		
	}			
	#--- OUTPUT
	res <- list(M_SD=M_SD, trafo_items=trafo_items, trafo_persons=trafo_persons,
					N_common=N_common)
	return(res)
}
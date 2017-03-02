
#***********************************************
# session info
Rsessinfo <- function(){
    si <- Sys.info()
    si2 <- utils::sessionInfo()
    paste0( si2$R.version$version.string , " " , si2$R.version$system 
             , " | nodename = " , si["nodename"] , " | login = " , si["login"] )
            }
#************************************************


#************************************************
# print CALL in summary										
print_CALL <- function(CALL){					
	cat("\n\nCall:\n", paste(deparse(CALL), sep = "\n", collapse = "\n"), 
				"\n\n", sep = "")	
							}
#************************************************							

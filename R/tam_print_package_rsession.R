## File Name: tam_print_package_rsession.R
## File Version: 0.03
## File Last Change: 2017-09-16 13:41:10

tam_print_package_rsession <- function(pack)
{
	for (pp in pack){
		cat( tam_packageinfo(pack=pp) , "\n" )	
	}
	cat( tam_rsessinfo() , "\n\n")	
}
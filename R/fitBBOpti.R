fitBBOpti <- function(filename=NULL){

	win <- if(.Platform$OS.type == "windows") TRUE else FALSE
	# if(!menus)win <- FALSE
	
	if(is.null(filename)){
		if(win){
			wd <- winDialog("okcancel", "Please select the CSV file with leaf gas exchange data.
			The BBOpti model will be fit using these data.\n
			Press OK to proceed.")
			
			if(wd == "CANCEL"){
				stop("Cancelled.\n")
			}
			# else...
			filename <- file.choose()
			
		} else {  # Mac/ unix.
			filename <- readline("Please enter name of the CSV file (without quotes): \n > ")
			if(!file.exists(filename)){
				message("File does not exist. These are CSV files in the current working directory:")
				print(list.files(pattern="\\.csv$", ignore.case=TRUE))
				stop()
			}
	
		}
	
	}
	
	dat <- try(read.csv(filename))
	# if(inherits(dat, "try-error"))stop("Please fix your dataset - it cannot be read using read.csv")

	A <- dat[,grep("^a$", tolower(names(dat)))]
	g <- dat[,grep("^gs$", tolower(names(dat)))]
	D <- dat[,grep("^d$", tolower(names(dat)))]
	if(length(D) == 0)D <- dat[,grep("^vpd$", tolower(names(dat)))]
	Ca <- dat[,grep("^ca$", tolower(names(dat)))]

	# message confirming all variables read. 
	if(length(A)>0 & length(g)>0 & length(D)>0 & length(Ca)>0){
		message("Successfully read: ", length(A), " rows of data.\n")
	} else {
		stop("One or more variables not read. Make sure D,gs,A and Ca are in the dataset.")
	}

	# Options: fit with or without a g0.
	# readline("Include g0 in the model fit (\'y\'), or assume g0 = 0.0 (\'n\')?")
	if(win){
		useg0 <- winDialog("yesno", "Include g0 in the model fit (\'Yes\'),\n or assume g0 = 0.0 (\'No\')?")
	} else {
		useg0 <- tolower(readline("Include g0 in the model fit (Y),\n or assume g0 = 0.0 (N)?\n (Y/N) > "))
		if(useg0 == "y")useg0 <- "YES"
		if(useg0 == "n")useg0 <- "NO"
	}

	# nls fit here...
	if(useg0 == "YES"){
		nlsfit <- nls(g ~ g0 + 1.6*(1+g1/sqrt(D))*(A/Ca), 	
			start=list(g0=0.02, g1=8))
		G0 <- coef(nlsfit)[[1]]
		G1 <- coef(nlsfit)[[2]]
		R2 <- 1-var(residuals(nlsfit))/var(g)
		ci <- suppressMessages(confint(nlsfit))
		G0ci <- c(ci[1],ci[3])
		G1ci <- c(ci[2],ci[4])
	}
	if(useg0 == "NO"){
		nlsfit <- nls(g ~ 1.6*(1+g1/sqrt(D))*(A/Ca), 	
			start=list(g1=8))
		G0 <- 0
		G1 <- coef(nlsfit)[[1]]
		R2 <- 1-var(residuals(nlsfit))/var(g)
		ci <- suppressMessages(confint(nlsfit))
		G0ci <- c(0,0)
		G1ci <- c(ci[1],ci[2])
	}

	if(win){
		savepred <- winDialog("yesno", "\nSave fitted values of the BBOpti model in a copy of the dataset?\n
		A CSV file will be written to the same folder as the original data.")
	} else {
		savepred <- readline("Save fitted values of the BBOpti model in a copy of the dataset?\n
		A CSV file will be written to the same folder as the original data.\n (Y/N) > ")
		if(savepred == "y")savepred <- "YES"
		if(savepred == "n")savepred <- "NO"
	}
	
	if(savepred == "YES"){
		gpred <- predict(nlsfit, data.frame(A=A,gs=g,Ca=Ca,D=D))
		dat$gBBOpti <- gpred
		fn <- paste(gsub("\\.csv$","",filename),"_R.csv")
		write.csv(dat, fn, row.names=FALSE)
	}
	
	# nice summary of coefficients and so on here...
	cat("\n")
	cat("----------------------------------------------\n")
	cat("Results:\n")
	cat("Fit of BBOpti model : gs = g0 + 1.6*(1 + g1/sqrt(D)) * (A/Ca).\n\n")
	if(useg0=="NO")cat("G0 assumed zero.\n")
	if(useg0=="YES"){
	cat("g0 =",round(G0,4),"  -- 95% CI =",round(G0ci[1],4),"-",round(G0ci[2],4),"\n")
	}
	cat("g1 =",round(G1,4),"  -- 95% CI =",round(G1ci[1],4),"-",round(G1ci[2],4),"\n\n")
	cat("R2 =", round(R2,3)," ,RMSE =", round(summary(nlsfit)$sigma,4),"\n")

	return(invisible(nlsfit))
}



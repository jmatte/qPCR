curve1 <- read.table("qPCRDataJPM.txt") # is possible to add name but getting the length of the table?

# the data need to be clean but after this get the log to find the linear portion of the curve
# remove any negative value, 

cyc <- length(curve1) - 1 #number of cycles minus the names (in this case only 1 that is the cell)
curves <- nrow(curve1)


log_Data <- matrix(nrow = curves, ncol = cyc)
    
for (cur in 1:curves){
    log_Data[cur,] <- rapply(curve1[cur,1:cyc+1],function(n){
        if(n<0){
            n<-10
        } else {
            n<-log10(n)
        }
    })
}

a <- 1
b <- cyc
minpoints <<- 4 # min number of points used in the curve

# number of possible lines between a and 5 with a min sep of minpoints
it <- ((b-minpoints+a)*(b-minpoints))/2
Cycles <- a:b # vector for the number of cycles.

# number of combination
# 50 cycles is 2+3+4+5...+49 itcal2<-((b1-3+a1)*(b1-3))/2+1 considering space of 3


Init_Value <- function (clean_Data){ #curv is a vector with the PCR data in LOG
    # matrix to record the result
    Result_Data <<- matrix(data = NA, nrow = curves, ncol = 6, 
                          dimnames = list(NULL, c("a","b","R^2","slope","C","Eff%")))
    for (cur in 1:curves){
        # matrix to record the R^2 of combinations
        bestline <- matrix(data = NA, nrow = it, ncol = 5)
        count <- 1
        for (i in a:b) {
            for (j in b:a) {
                if ( (i < j) & (j-i >= minpoints) ) {
                    Datatest <- as.data.frame(t(rbind(Cycles[(i):(j)],clean_Data[cur,(i):(j)])))
                    bestline[count,1:2] <- (c(i,j))
                    Temp_Data_lm <- lm(Datatest[,1] ~ Datatest[,2], 
                                       data = Datatest)
                    bestline[count,3] <- summary(Temp_Data_lm)$r.squared
                    bestline[count,4] <- Temp_Data_lm$coefficients[[2]]
                    bestline[count,5] <- Temp_Data_lm$coefficients[[1]]
                    count <- count + 1
                } 
            }
        }
        Result_Data[cur,1:5] <<- bestline[which.max(bestline[,3]),]
        Result_Data[cur,6] <<- 10^(1/Result_Data[cur,4])-1
    }
}


Init_Value(log_Data)

head(Result_Data)
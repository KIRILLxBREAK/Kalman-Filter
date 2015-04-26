sigma <- 1
sigma0 <- 1
T <- 1

F1 <- c(1, T)
F2 <- c(0, 1)
F <- rbind(F1, F2)

P1 <- c(sigma*sigma, sigma*sigma/T)
P2 <- c(sigma*sigma/T, 2*sigma*sigma/(T*T))
P <- rbind(P1, P2)

H <- c(1, 0)
R <- sigma
Q <- matrix(c(0,0,0,0), 2, 2)
q = 0
I = matrix(c(1, 0, 0, 1), 2, 2)

z <- vector("numeric", length = 100)
x <- vector("numeric", length = 100)

z[,1] <- c(1, 8000)
z[,2] <- c(2, 1)

for (i in 1:100) {
    if (i==1) {
        x[,1] <- z[,1]

    } else if (i==2) {
        x[1,2] <- z[1,2]
        x[2,2] <- (x[1,2]-x[1,1])/T

    } else

        ##prediction
        z[,i] <- F*z[,i-1]
        x[,i] <- F*x[,i-1]
               
        P <- F*P*F
               
        
        ##correction
        K <- P*H/(H*P*H+R)
          
        x[,i] <- x[,i] + K*(z[1,i] - H*x[,i])
        P <- (I-K*H)*P
    }   
}


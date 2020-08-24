# Compare different functions for homeostatic activation

setwd("~/Desktop/HIV_model")

source("grind.R")
colours <- c("blue", "darkgreen", "red", "darkorange","darkmagenta") #,"gold","darkorchid","aquamarine","deeppink","gray")

############## The loop model - Convex HA ###########################################
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    
    kt = k*exp(-alpha*t)
    Ket = Ke + eps_Ke*t
    wit = wi + eps_wi*t
    w = max(w0*(1-(R+Tc+I)/(K+Tc+R+I)), 0)
    wa = wit*I/(ha + I)
    beta_i = beta*(I/(Ki+I))
    dI = 2*(wa+w)*beta_i*R - di*I - kt*I*E
    dR = r*Tc - (wa+w)*R - dr*R
    dTc = 2*(wa+w)*(1-beta_i)*R - r*Tc - dt*Tc
    dE = p*(E*I)/(Ket+E+I) - de*E
    
    return(list(c(dI, dR, dTc, dE)))  
  }) 
}

p <- c(r=1/5, w0=1/150, dt=1/150, dr=1/150, d=1/150, K=1e6, di=0.1, 
       beta=0, Ki=1e3, p=1.1, k=1e-4, de=0.1, Ke=100, wi=1/150,
       ha=3.3e4, alpha=0, eps_Ke=0, eps_wi=0)

p["K"] <- 1/3*1e6
p["w0"] <- 4/150
p["r"] <- 33/50 #19701/15000  # 2*99*1/150 # #99*(1/150) + 99*p["dr"]  # w_ss - dt0 (1/150 = wrt in ss)
p["ha"] <- 3.3e4

s <- c(I=0, R=1e5, Tc=1e2, E=0)
s <- run(2000, log="y", ymin = 1)
s <- newton(s)

# Initiate infection and immune response
s["I"] <-1
s["E"] <- 1
p["beta"] <- 1
p["Ki"] <-3.2e3

# Early infection
sol100 <- run(100, log="y", table=T, ymin=1,
              tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                        "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"))

abline(h=3, lty=2)

# Look at viral replication rate during initial increase
calc_rho <- sol100[sol100$I<=max(sol100$I[1:40])& sol100$time<=sol100$time[sol100$I==max(sol100$I[1:40])],]
calc_rho$I <- log10(calc_rho$I)

plot(calc_rho$time, calc_rho$I, type = "p", ylab = "I cell density", xlab = "time")
abline(0,1.5, col="red", lty=2)

(log10(sol100$I[3])-log10(sol100$I[1]))/2

######### Ten years infection
p["ha"] <- 3.3e4
p["Ki"] <- 3.2e3
v_sol_HA <- run(3650, log="y",table = T, tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                                                 "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"),
              main="Homeostatic activation")
abline(h=3, lty=2)
abline(h=200e3, lty=3)

p["ha"] <- 1e3
p["Ki"] <- 3.2e3
v_sol_IA <- run(3650, log="y",table = T, tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                                                 "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"),
              main="Immune activation")
abline(h=3, lty=2)
abline(h=200e3, lty=3)

############################## Linear HA  ###########################################
model <- function(t, state, parms) {
  
  with(as.list(c(state,parms)), {
    
    kt = k*exp(-alpha*t)
    Ket = Ke + eps_Ke*t
    wit = wi + eps_wi*t
    w = max(w0*(1-(R+Tc+I)/K),0)
    wa = wit*I/(ha + I)
    beta_i = beta*(I/(Ki+I))
    dI = 2*(wa+w)*beta_i*R - di*I - kt*I*E
    dR = r*Tc - (wa+w)*R - d*R
    dTc = 2*(wa+w)*(1-beta_i)*R - r*Tc - d*Tc
    dE = p*(E*I)/(Ket+E+I) - de*E
    
    return(list(c(dI, dR, dTc, dE)))  
  }) 
}

p <- c(r=1/5, w0=1/150, d=1/150, d=1/150, K=1e6, di=0.1, 
       beta=0, Ki=1e3, p=1.1, k=1e-4, de=0.1, Ke=100, wi=1/150,
       ha=3.3e4, alpha=0, eps_Ke=0, eps_wi=0)

p["K"] <- 4/3*1e6
p["w0"] <- 4/150
p["r"] <- 33/50 #19701/15000 # #99*(1/150) + 99*p["dr"]  # w_ss - dt0 (1/150 = wrt in ss)
p["ha"] <- 3.3e4

s <- c(I=0, R=1e6, Tc=1e2, E=0)
s <- run(2000, log="y", ymin = 1)
s <- newton(s)

s["E"] <- 1
s["I"] <- 1
p["beta"] <- 1
p["Ki"] <- 3.25e3

# Early infection
sol100 <- run(100, log="y", table=T, tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                                               "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"))

abline(h=3, lty=2)

# Look at viral replication rate during initial increase
calc_rho <- sol100[sol100$I<=max(sol100$I)& sol100$time<=sol100$time[sol100$I==max(sol100$I)],]
calc_rho$I <- log10(calc_rho$I)

plot(calc_rho$time, calc_rho$I, type = "p", ylab = "I cell density", xlab = "time")
abline(0,1.5, col="red", lty=2)

(log10(sol100$I[3])-log10(sol100$I[1]))/2

######### Ten years infection
p["ha"] <- 3.3e4
p["Ki"] <- 3.25e3
l_sol_HA <- run(3650, log="y", table=T, ymin=1, tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                                                     "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"))
abline(h=3, lty=2)
tail(sol)

p["ha"] <- 1e3
p["Ki"] <- 3.25e3

l_sol_IA <- run(3650, log="y", table=T, ymin=1, tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                                                          "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"))
abline(h=3, lty=2)

################## Concave HA ########################################
# The loop model
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    
    kt = k*exp(-alpha*t)
    Ket = Ke + eps_Ke*t
    wit = wi + eps_wi*t
    w = max(w0*(1/(1+((Tc+R+I)/K)^n)), 0)
    wa = wit*I/(ha + I)
    beta_i = beta*(I/(Ki+I))
    dI = 2*(wa+w)*beta_i*R - di*I - kt*I*E
    dR = r*Tc - (wa+w)*R - d*R
    dTc = 2*(wa+w)*(1-beta_i)*R - r*Tc - d*Tc
    dE = p*(E*I)/(Ket+E+I) - de*E
    
    return(list(c(dI, dR, dTc, dE)))  
  }) 
}

p <- c(r=1/5, w0=1/150, d=1/150, d=1/150, K=1e6, di=0.1, 
       beta=0, Ki=1e3, p=1.1, k=1e-4, de=0.1, Ke=100, wi=1/150,
       ha=3e4, alpha=0, eps_Ke=0, eps_wi=0)

p["n"] <- 2
p["K"] <- 1/sqrt(3)*1e6
p["w0"] <- 4/150
p["r"] <- 33/50 #19701/15000  # 2*99*1/150 # #99*(1/150) + 99*p["dr"]  # w_ss - dt0 (1/150 = wrt in ss)
p["ha"] <- 3.3e4

s <- c(I=0, R=1e6, Tc=1e2, E=0)
s <- run(2000, log="y", ymin = 1)
s['Tc']/sum(s)*100
s <- newton(s)
sum(s)

# Initiate infection and immune response
s["I"] <-1
s["E"] <- 1
p["beta"] <- 1
p["Ki"] <-3.25e3

# Early infection
sol100 <- run(100, log="y", table=T, ymin=1,
              tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                        "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"))

abline(h=3, lty=2)

# Look at viral replication rate during initial increase
calc_rho <- sol100[sol100$I<=max(sol100$I[1:40])& sol100$time<=sol100$time[sol100$I==max(sol100$I[1:40])],]
calc_rho$I <- log10(calc_rho$I)

plot(calc_rho$time, calc_rho$I, type = "p")
abline(0,1.5, col="red", lty=2)

(log10(sol100$I[3])-log10(sol100$I[1]))/2

# Ten years infection
p["ha"] <- 3.3e4
p["Ki"] <- 3.25e3
c_sol_HA <- run(3650, log="y",table = T, tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                                              "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"))
abline(h=3, lty=2)

p["ha"] <- 1e3
p["Ki"] <- 3.25e3
c_sol_IA <- run(3650, log="y",table = T, tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                                                   "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"))
abline(h=3, lty=2)


################# Plot figure #########################################
plot(v_sol_HA$time, v_sol_HA$R, type="l", col="blue", lwd=2,
     log="y", xlim=c(0,1500), ylim=c(1.8e5, 4e5),
     main="Comparing different homeostatic activation functions", xlab="Time", ylab="R cell density (cells/mL)", 
     cex.axis=1.2, cex.lab=1.2, cex.main=1.2)
lines(l_sol_HA$time, l_sol_HA$R, col="red", lwd=2)
#lines(l_sol_IA$time, l_sol_IA$R, col="red", lwd=2, lty=2)
lines(c_sol_HA$time, c_sol_HA$R, col="darkgreen", lwd=2)
#lines(c_sol_IA$time, c_sol_IA$R, col="darkgreen", lwd=2, lty=2)
abline(h=200e3, lty=3, lwd=2)
legend("topright", 
       legend = c("Convex HA", "Linear HA", "Concave HA"),
       col=c("blue", "red", "darkgreen"),
       lty=1,
       lwd=2,
       cex=0.9)

# This scripts models HIV infection. The within host model of van Dorp et al. 2020 is 
# simplified to one viral strain and one immmune response. It is subsequently adapted
# to attempt to model the chronic phase of HIV infection. Specifically focussing on
# linking infection explicitly to cell division.

# Created by: Alice Driessen
# Date: 19-05-2020

setwd("~/Desktop/HIV_model")

source("grind.R")
colours <- c("blue", "darkgreen", "red", "darkorange","darkmagenta") #,"gold","darkorchid","aquamarine","deeppink","gray")

# The loop model
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
s['Tc']/sum(s)*100
s <- newton(s)
sum(s)

# Depletion of cells in healthy state as with chemotherapy of stem cell transplantation
sol_Depl <- run(2000, table=T, after = "if(t==500){state['R'] <- 100000}", log="y", ymin=1e3)

plot(sol_Depl$time, sol_Depl$R, type="l", col="blue", lwd=2, log="y", ylim=c(1e3,1e6),
     ylab = "Cell Density (cells/mL)", xlab="Time (days)", main="Depletion in uninfected model",
     cex.axis=0.9, cex.lab=0.9, cex.main=1.2)
lines(sol_Depl$time, sol_Depl$Tc, lwd=2, col="darkgreen")
legend("topright", legend = c("R", "T"), lty=1, lwd=2, col=colours, cex=0.9)

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

plot(calc_rho$time, calc_rho$I, type = "p", ylab = "I cell density", xlab = "Time (days)")
abline(0,1.5, col="red", lty=2)

(log10(sol100$I[3])-log10(sol100$I[1]))/2 # Should be 1-1.5

######### Ten years infection ####################################
p["ha"] <- 3.3e4
p["Ki"] <- 3.2e3
sol_HA <- run(3650, log="y",table = T, tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                                              "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"),
           main="Homeostatic activation")
abline(h=3, lty=2)
abline(h=200e3, lty=3)
tail(sol)

plot(sol_HA$time, sol_HA$R, log="y", type = "l", lwd=2, col="blue", xlim=c(0,150),
     ylab = "Cell Density (cells/mL)", xlab="Time (days)", main="Homeostatic activation", ylim = c(1, 1e6),
     cex.axis=1.2, cex.lab=1.2, cex.main=1.2)
lines(sol_HA$time, sol_HA$Tc, col="darkgreen", lwd=2)
lines(sol_HA$time, sol_HA$I, col="red", lwd=2)
lines(sol_HA$time, sol_HA$E, col="darkorange", lwd=2)
lines(sol_HA$time, sol_HA$`%Tc+I`, col="darkmagenta", lwd=2, lty=2)
abline(h=3, lty=1, col="darkgrey")
abline(h=200e3, lty=3, lwd=2)
legend("topright", legend = c("R", "T", "I", "E", "%T+I"), lty=c(1,1,1,1,2), lwd=2, col=colours[1:5], cex=0.9)
mtext(text="C", side=3, line=2, at=-20, cex=2.5)

plot(x=sol_HA$time, y=sol_HA$R/1000, col="blue", type="l", 
     xlab = "Time (days)", ylab = "R cells/microliter", main = "Homeostatic activation",
     lwd=2, xlim=c(0,1500), cex.axis=1.2, cex.lab=1.2, cex.main=1.2)
abline(h=200, lty=2)

sum_HA <- sum(sol_HA[3651,c("R", "Tc", "I")])
turnover_HA <- p["d"]*sum(sol_HA[3651,c("R", "Tc")])/sum_HA + (p["k"]*sol_HA[3651,"E"]+p["di"])*sol_HA[3651,"I"]/sum_HA

p["ha"] <- 1e3
p["Ki"] <- 3.2e3
sol_IA <- run(3650, log="y",table = T, tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                                                 "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"),
              main="Immune activation")
abline(h=3, lty=2)
abline(h=200e3, lty=3)

sum_IA <- sum(sol_IA[3651,c("R", "Tc", "I")])
turnover_IA <- p["d"]*sum(sol_IA[3651,c("R", "Tc")])/sum_IA + (p["k"]*sol_IA[3651,"E"]+p["di"])*sol_IA[3651,"I"]/sum_IA


plot(sol_IA$time, sol_IA$R, log="y", type = "l", lwd=2, col="blue", xlim=c(0,150),
     ylab = "Cell Density (cells/mL)", xlab="Time (days)", main="Homeostatic and immune activation", ylim = c(1, 1e6),
     cex.axis=1.2, cex.lab=1.2, cex.main=1.2)
lines(sol_IA$time, sol_IA$Tc, col="darkgreen", lwd=2)
lines(sol_IA$time, sol_IA$I, col="red", lwd=2)
lines(sol_IA$time, sol_IA$E, col="darkorange", lwd=2)
lines(sol_IA$time, sol_IA$`%Tc+I`, col="darkmagenta", lwd=2, lty=2)
abline(h=3, lty=1, col="darkgrey")
abline(h=200e3, lty=3, lwd=2)
legend("topright", legend = c("R", "T", "I", "E", "%T+I"), lty=c(1,1,1,1,2), lwd=2, col=colours, cex=0.9)
mtext(text="D", side=3, line=2, at=-20, cex=2.5)

p["ha"] <- 1e3
p["Ki"] <- 3.5e3
sol_IA_Ki <- run(3650, log="y",table = T, tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                                                    "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"),
                 main="Immune activation")

p["ha"] <- 1e3
p["Ki"] <- 2.3e3
sol_IA_KiLow <- run(3650, log="y",table = T, tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                                                    "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"),
                 main="Immune activation")

p["ha"] <- 3.3e4
p["Ki"] <- 2.3e3
sol_HA_KiLow <- run(3650, log="y",table = T, tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                                                       "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"),
                    main="Immune activation")


# R cells in cell/microL
plot(x=sol_HA$time, y=sol_HA$R/1000, col="blue", type="l", 
     xlab = "Time (days)", ylab = "R cells/microliter", main = "Comparing homeostatic and immune activation",
     lwd=2, xlim=c(0,1500), ylim = c(100,550), cex.axis=1.2, cex.lab=1.2, cex.main=1.2, log="y")
abline(h=200, lty=3, lwd=2)
lines(sol_IA$time, sol_IA$R/1000, col="red", lwd=2)
lines(sol_IA_Ki$time, sol_IA_Ki$R/1000, col="purple", lwd=2)
lines(sol_IA_KiLow$time, sol_IA_KiLow$R/1000, col="black", lwd=2)
lines(sol_HA_KiLow$time, sol_HA_KiLow$R/1000, col="green", lwd=2)
legend("topright",
       legend = c("Homeostatic activation, Ki=3200; rho=1.6", "Homeostatic and immune activation, Ki=3200; rho=1.8", 
                  "Homeostatic and immune activation, Ki=3500; rho=1.6",  "Homeostatic activation, Ki=2300; rho=1.9",
                  "Homeostatic and immune activation, Ki=2300; rho=2.1"), 
       col=c("blue", "red", "purple", "green", "black"), lty=1, cex=0.9, lwd=2)
mtext(text="A", side=3, line=2, at=-200, cex=2.5)


# Percentages of T + I cells in comparing immune activation and homeostatic activation
plot(x=sol_HA$time, y=(sol_HA$I+sol_HA$Tc)/(sol_HA$I+sol_HA$Tc+sol_HA$R)*100, col="blue", type="l", lty=2,
     xlab = "Time (days)", ylab = "% recently divided cells", main = "Comparing homeostatic and immune activation",
     lwd=2, xlim=c(0,1500), cex.axis=1.2, cex.lab=1.2, cex.main=1.2, ylim=c(1,11))
abline(h=3, lty=1)
lines(sol_IA$time, (sol_IA$I+sol_IA$Tc)/(sol_IA$Tc+sol_IA$R+sol_IA$I)*100, col="red", lwd=2, lty=2)
lines(sol_IA_Ki$time, (sol_IA_Ki$Tc+sol_IA_Ki$I)/(sol_IA_Ki$Tc+sol_IA_Ki$I+sol_IA_Ki$R)*100, col="purple", lwd=2, lty=2)
lines(sol_IA_KiLow$time, (sol_IA_KiLow$Tc+sol_IA_KiLow$I)/(sol_IA_KiLow$Tc+sol_IA_KiLow$I+sol_IA_KiLow$R)*100, col="black", lwd=2, lty=2)
lines(sol_HA_KiLow$time, (sol_HA_KiLow$Tc+sol_HA_KiLow$I)/(sol_HA_KiLow$Tc+sol_HA_KiLow$I+sol_HA_KiLow$R)*100, col="green", lwd=2, lty=2)
legend("topright",
       legend = c("Homeostatic activation, Ki=3200; rho=1.6", "Homeostatic and immune activation, Ki=3200; rho=1.8", 
                  "Homeostatic and immune activation, Ki=3500; rho=1.6",  "Homeostatic activation, Ki=2300; rho=1.9",
                  "Homeostatic and immune activation, Ki=2300; rho=2.1"), 
       col=c("blue", "red", "purple", "green", "black"), lty=2, cex=0.9, lwd=2)
mtext(text="B", side=3, line=2, at=-200, cex=2.5)

#abline(h=0.5*sol_microL$R[1], lty=2)
sol_microL[sol_microL$R<201&sol_microL$R>199,]

sol[sol$R<(0.5*sol$R[1]+1000) & sol$R>(0.5*sol$R[1]-1000),]

plot(sol$time, sol$Tc+sol$R, type = "l", xlim=c(0,100))

######### Include disease progession ##########################
p["Ki"] <- 5.5e3
p["ha"] <- 1e3

p["alpha"] <-  0 #3.5e-4
p["eps_Ke"] <- 0# 6
p["eps_wi"] <- 3.6e-5 #3.6e-5


solDP <- run(4500, log="y", table = T, ymin=1, 
             tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                       "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"))
abline(h=200e3)

plot(solDP$time, solDP$R, log="y", type = "l", lwd=2, col="blue",
     ylab = "Cell Density (cells/mL)", xlab="Time (days)", main="Increasing immune activation", ylim = c(1, 1e6),
     cex.axis=1.2, cex.lab=1.2, cex.main=1.2)
lines(solDP$time, solDP$Tc, col="darkgreen", lwd=2)
lines(solDP$time, solDP$I, col="red", lwd=2)
lines(solDP$time, solDP$E, col="darkorange", lwd=2)
lines(solDP$time, solDP$`%Tc+I`, col="darkmagenta", lwd=2, lty=2)
abline(h=3, lty=1, col="darkgrey")
abline(h=200e3, lty=3, lwd=2)
legend("topright", legend = c("R", "T", "I", "E", "%T+I"), lty=c(1,1,1,1,2), lwd=2, col=colours, cex=0.9)
mtext(text="A", side=3, line=2, at=-600, cex=2.5)

# R cells on linear scale
plot(solDP$time, solDP$R, type="l", col="blue",lwd=2, 
     ylab = "R Cell Density (cells/mL)", xlab="Time (days)", main="Increasing immune activation",
     cex.axis=1.2, cex.lab=1.2, cex.main=1.2)
abline(h=200e3, lty=2)
mtext(text="B", side=3, line=2, at=-600, cex=2.5)

solDP[solDP$R<200e3,][1,]
tail(solDP)

# First three months
plot(solDP$time, solDP$R, type="l", col="blue",lwd=2, xlim = c(0,120))
lines(solDP$time, solDP$Tc, col="green", lwd=2)
lines(solDP$time, solDP$I, col="red", lwd=2)

#################### Lymph and Blood compartments ##########################################################
p["alpha"] <-  0 #3.5e-4
p["eps_Ke"] <- 0# 6
p["eps_wi"] <- 0 #3.6e-5

p["ha"] <- 3.3e4
p["Ki"] <- 3.2e3

sol_HAART <- run(1500, log="y", ymin=1, table=T, after = "if(t==500){parms['beta']<- 0}",
                 tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                           "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"))

plot(sol_HAART$time, sol_HAART$R, log="y", type = "l", lwd=2, col="blue",
     ylab = "Cell Density (cells/mL)", xlab="Time (days)", main="HAART", ylim = c(0.2, 1e6),
     cex.axis=1.2, cex.lab=1.2, cex.main=1.2)
lines(sol_HAART$time, sol_HAART$Tc, col="darkgreen", lwd=2)
lines(sol_HAART$time, sol_HAART$I, col="red", lwd=2)
lines(sol_HAART$time, sol_HAART$E, col="darkorange", lwd=2)
lines(sol_HAART$time, sol_HAART$`%Tc`, col="darkmagenta", lwd=2, lty=2)
abline(h=2, lty=4, col="darkgrey")
legend("topright", legend = c("R", "T", "I", "E", "%T"), lty=c(1,1,1,1,2), lwd=2, col=colours, cex=0.9)
mtext(text="A", side=3, line=2, at=-200, cex=2.5)


sol_LB <- run(1500, log="y", ymin = 1, table=T,
           tweak = c("nsol$fr = 0.02 - (0.01*nsol$I)/(1e3+nsol$I)",
                     "nsol$ft = 0.02 + (0.09*nsol$I)/(1e3+nsol$I)",
                     "nsol$Rb = nsol$fr*nsol$R",
                     "nsol$Rl = (1-nsol$fr)*nsol$R",
                     "nsol$Tb = nsol$ft*nsol$Tc",
                     "nsol$Tl = (1-nsol$ft)*nsol$Tc",
                     "nsol$`%T_b` = nsol$Tb/(nsol$Tb+nsol$Rb)*100",
                     "nsol$Tot_b = nsol$Tb+nsol$Rb",
                     "nsol$`%B` = nsol$Tot_b/(nsol$Tc+nsol$R)*100",
                     "nsol$`%T_l` = nsol$Tl/(nsol$Tl+nsol$R)*100"),
           after = "if(t==500){parms['beta']<- 0}")
abline(h=2, lty=2)

plot(sol_LB$time, sol_LB$Rb, type="l", col="blue", lwd="2", log="y", ylim = c(0.2, 1e6),
     ylab = "Cell Density (cells/mL)", xlab="Time (days)", main="HAART - blood", 
     cex.axis=1.2, cex.lab=1.2, cex.main=1.2)
lines(sol_LB$time, sol_LB$Tb, lwd=2, col="darkgreen")
lines(sol_LB$time, sol_LB$`%T_b`, lwd=2, col="darkmagenta", lty=2)
lines(sol_LB$time, sol_LB$Tot_b, lwd=2, col="red")
lines(sol_LB$time, sol_LB$`%B`, lty=2, lwd=2, col="gold")
abline(h=sol_LB[1,"Tb"], lty=4, lwd=1, col="darkgrey")
abline(h=sol_LB[1,"%T_b"], lty=4, lwd=1, col="darkgrey")
legend("topright", legend=c("Rb", "Tb", "Tot_b", "%Tb", "%Tot_b"),
       col=c("blue", "darkgreen", "red", "darkmagenta", "gold"), lty=c(1,1,1,2,2),
       cex=0.9, lwd=2)
mtext(text="B", side=3, line=2, at=-200, cex=2.5)

plot(sol_LB$time, sol_LB$`%B`, lwd=2, col="red", lty=3, type='l', ylim = c(0,30))
abline(h=2, lty=4)

plot(sol_LB$time, sol_LB$Rb, col="blue", type="l", lty=2, xlim = c(0,100))
abline(h=sol_LB[1,"Rb"], lty=3)

# Lymph compartment
colnames(sol_LB)
sol_LB$`%T_l` <- apply(sol_LB, 1, function(row) {row['Tl']/(row['Tl']+row['R'])*100})

plot(sol_LB$time, sol_LB$Rl, type="l", col="blue", lwd="2", log="y", 
     ylim=c(1e-1,1e6), xlab = "Time (days)", ylab = 'density', main="Lymphoid tissue")
lines(sol_LB$time, sol_LB$Tl, lwd=2, col="green")
lines(sol_LB$time, sol_LB$`%T_l`, lwd=2, lty=2, col='aquamarine')
abline(h=1, lty=2, col='grey')
abline(h=sol_LB[1,"Tl"], lty=2)
legend("topright", legend=c("Rl", "Tl", "%T_l"),
       col=c("blue", "green", "black"), lty=c(1,1,1),
       cex=0.9)

head(sol_LB)

#################### Higher death rate -- Grossman 2002 #################################
p["alpha"] <-  0 #3.5e-4
p["eps_Ke"] <- 0# 6
p["eps_wi"] <- 0 #3.6e-5

p["K"] <- 1/3*1e6
p["w0"] <- 4/150
p["wi"] <- 1/150

p["dt"] <- 49/4350  #= 0.0112644
p["r"] <- 2891/4350 #= 0.664598

s <- c(I=0, R=1e5, Tc=1e2, E=0)
s <- run(2000, log="y", ymin = 1)
s['Tc']/sum(s)*100
s <- newton(s)
sum(s)

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

plot(calc_rho$time, calc_rho$I, type = "p", ylab = "I cell density", xlab = "Time (days)")
abline(0,1.5, col="red", lty=2)

(log10(sol100$I[3])-log10(sol100$I[1]))/2

# Ten years infection 
p["ha"] <- 3.3e4
p["Ki"] <- 3.2e3
sol_HA_dt <- run(3650, log="y",table = T, tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                                                 "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"),
              main="Homeostatic activation")
abline(h=3, lty=2)
abline(h=200e3, lty=3)
tail(sol_HA_dt)
head(sol_HA_dt[sol_HA_dt$R<2e5,])

p["ha"] <- 1e3
sol_IA_dt <- run(3650, log="y",table = T, tweak = c("nsol$`%Tc` = nsol$Tc/(nsol$Tc+nsol$R)*100",
                                                 "nsol$`%Tc+I` = (nsol$Tc+nsol$I)/(nsol$Tc+nsol$R+nsol$I)*100"),
              main="Homeostatic activation")
abline(h=3, lty=2)
abline(h=200e3, lty=3)
tail(sol_IA_dt)
head(sol_IA[sol_IA$R<2e5,])


plot(sol_IA_dt$time, sol_IA_dt$R, type="l", log="y", col="red", lty=2, lwd=2, 
     xlim=c(0,1500), ylim=c(1.9e5, 4e5),
     main="Higher death rate for recently activated cells", xlab="Time (days)", ylab="R cell density (cells/mL)", 
     cex.axis=1.2, cex.lab=1.2, cex.main=1.2)
lines(sol_IA$time, sol_IA$R, col="blue", lty=2, lwd=2)
lines(sol_HA$time, sol_HA$R, col="blue", lwd=2)
lines(sol_HA_dt$time, sol_HA_dt$R, col="red", lwd=2)
abline(h=2e5, lty=3, lwd=2)
legend("topright", legend=c("Equal death rate", "Increased T death rate", "Equal death rate + IA", "Increased T death rate + IA"), 
      col=c("blue", "red", "blue", "red"), lty=c(1,1,2,2), lwd=2, cex=0.9)


################## Activation rates, infection chance, loss and gain #################################
sol$w <- p["w0"]*(1-(sol$I+sol$R+sol$Tc)/(sol$I+sol$R+sol$Tc+p["K"]))
sol$wit <- p["wi"] + p["eps_wi"]*sol$t
sol$wa <- sol$wit*sol$I/(sol$I+p["ha"])
sol$betaI <- p["beta"]*sol$I/(p["Ki"]+sol$I)

plot(sol$time, sol$w, type = "l")
plot(sol$time, sol$wa, type = "l")
plot(sol$time, sol$wa+sol$w, type = "l")
abline(h=1/150, lty=2, col="red")
tail(sol)
tail(sol$w+sol$wa)

plot(sol$time, sol$betaI, type="l")
lines(sol$time, 0.5-p["d"]*(sol$R+sol$Tc)/(2*sol$R*(sol$wa+sol$w)), lty=2, col="red")
abline(h=0.5-p["d"]/(2*p['r']), lty=2, col="green")

plot(sol$time, 2*(1-sol$betaI)*p["r"]/(p["r"]+p["d"])-1, type = "l", ylim=c(-1,2))
lines(sol$time, 1-(2*(1-sol$betaI)*p["r"]/(p["r"]+p["d"])-1), lty=2)
2*(1-sol$betaI[3000])*p["r"]/(p["r"]+p["d"])-1

# Total population loss vs gain
plot(sol$time,  sol$R*(sol$wa+sol$w)*(1-2*sol$betaI), type="l")
lines(sol$time, p["d"]*(sol$R+sol$Tc), lty=2)

plot(sol$time, (sol$R*(sol$wa+sol$w)*(1-2*sol$betaI))-(p["d"]*(sol$R+sol$Tc)), 
     type="l", ylim=c(-5000, 100))
abline(h=0, lty=2)

# R population loss vs gain
plot(sol$time,  2*sol$R*(sol$wa+sol$w)*(1-sol$betaI)*p["r"]/(p["r"]+p["d"]), type="l")
lines(sol$time, (sol$w+sol$wa+p["d"])*(sol$R), lty=2)

plot(sol$time, (2*sol$R*(sol$wa+sol$w)*(1-sol$betaI)*p["r"]/(p["r"]+p["d"]))-((sol$w+sol$wa+p["d"])*(sol$R)), 
     type="l", ylim=c(-5000, 100))
abline(h=0, lty=2)

# T population loss vs gain
plot(sol$time,  2*sol$R*(sol$wa+sol$w)*(1-sol$betaI), type="l")
lines(sol$time, (p["r"]+p["d"])*(sol$Tc), lty=2)

# Plot scaling of Wrt
test <- data.frame("tot" = seq(1e4, 1e7, 1000))
test$K1 <- p["K"]
test$f1 <- apply(test, 1, function(row){1-row["tot"]/(row["tot"]+row["K1"])})
test$w01 <- p["w0"]
test$w1 <- test$f1*unique(test$w01)

plot(test$tot, test$w1, type="l", ylim=c(0,unique(test$w01)), xlim = c(1e4, 1e7),
     xlab= "Total cell density", ylab="homeostatic activation", main = "Convex scaling")
abline(v=p["K"], lty=2)
abline(h=0.5*p["w0"], lty=2)
abline(h=1/150, col="red")
abline(v=sum(s), col="red")
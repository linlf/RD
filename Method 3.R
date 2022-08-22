# Input : Baseline risk and CI (BR BRLL and BRUL), relative risk and CI (RR, RRLL, RRUL)

BR <- 0.053; BRLL <- 0.039; BRUL <- 0.072
RR <- 0.77; RRLL <- 0.63; RRUL <- 0.94
N <- 10000 # No. of simulations
# Get shape parameters of beta dist 
BR_variance <- ((BRUL - BRLL)/3.92)^2
alpha <- ((1 - BR) / BR_variance - 1 / BR) * BR ^ 2
beta <- alpha * (1 / BR - 1)
# Simulation
sim_BR <- rbeta(N, alpha, beta, ncp = 0)
sim_RR <- rlnorm(N, meanlog = log(RR), sdlog = (log(RRUL) - log(RRLL))/3.92)
RD <- sim_BR*(sim_RR - 1)*1000 # Multiplying by 1,000 scales RD to per 1,000 patients
# Output
mean(RD); median(RD)
round(quantile(RD, c(0.025, 0.975)), 2)
hist(RD)

# Input : Baseline risk available as a range (BRLL and BRUL), relative risk and CI (RR, RRLL, RRUL)
BRLL <- 0.039; BRUL <- 0.072
RR <- 0.77; RRLL <- 0.63; RRUL <- 0.94
N <- 10000 # No. of simulations

sim_BR <- runif(N, min = BRLL, max = BRUL)
sim_RR <- rlnorm(N, meanlog = log(RR), sdlog = (log(RRUL) - log(RRLL))/3.92)
RD <- sim_BR*(sim_RR - 1)*1000 # 1,000 scales RD to per 1,000 patients

mean(RD); median(RD)
round(quantile(RD, c(0.025, 0.975)), 2)
tiff(filename = "RD Histogram2.tif", width = 2800, height = 1000,
  pointsize = 10, res = 300)
hist(RD)
dev.off()
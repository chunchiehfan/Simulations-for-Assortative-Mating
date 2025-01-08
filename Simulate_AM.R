
# This is a script to simulate scenarios about AM

# package for tetrachortic correlations
require(psych)


# cohort size
N = 10000 

# age strata
ageset <- c(20,40,60,80,Inf)


# Scenario 1: 
# No censoring
# There are prevalence difference across age strata
#     prevalence scaled by age monotonically with the scalar slope b

assortvec <- c(0.1,0.2,0.3,0.4,0.5)
slopevec <- c(5,10,15,20) # difference in the first order
df <- data.frame()
for (assort in assortvec){
for (s in slopevec){
for (iiter in 1:100){
print(paste0(iiter))
agevec.1 <- runif(N, min = 25, max=90)
agevec.2 <- agevec.1 + runif(N, min = -5, max=5) # random assign the age gap, at most 10
prev.1 <- 0.05 + as.matrix(poly(agevec.1, 2))%*%matrix(c(s,-2),2,1) # assuming weibull with lifetime risk at 0.14
prev.2 <- 0.05 + as.matrix(poly(agevec.2, 2))%*%matrix(c(s,-2),2,1)
liab.1 <- rnorm(N)
liab.2 <- (assort)*liab.1 + sqrt(1 - assort**2)*rnorm(N)
casevec.1 <- liab.1 > qnorm(prev.1, lower=F)
casevec.2 <- liab.2  > qnorm(prev.2, lower=F)
rvec0 <- c()
rvec1 <- c()
prev0 <- c()
for (ibin in 2:length(ageset)){
  inset <- (agevec.1 >= ageset[ibin-1]) & (agevec.1 < ageset[ibin])
  prev0 <- c(prev0, mean(c(prev.1[inset],prev.2[inset])))
  rvec0 <- c(rvec0, cor(liab.1[inset],liab.2[inset]))
  rvec1 <- c(rvec1, tetrachoric(matrix(table(casevec.1[inset],casevec.2[inset]),2,2))$rho)
}
tmp <- data.frame(liab=rvec0, tetra=rvec1, prev=prev0, iter=iiter, assort=assort, age1 = ageset[1:(length(ageset)-1)], age2 = ageset[2:length(ageset)])
df <- rbind(df, tmp)
}}} 

df$ageset <- paste0(df$age1, ' ~ ', df$age2)
ggplot(df, aes(prev, tetra, fill=factor(assort))) + 
	geom_point(shape=21,alpha=0.5,size=3,colour='grey') + 
	geom_hline(yintercept=assortvec,linetype=2) + 
	theme_bw() + 
	labs(x = 'prevalence', y='tetrachoric correlations', fill='Level of Assortment')
ggsave('Simulations_birth_year_dependent_prev.pdf')

ggplot(subset(df, assort==0.5), aes(ageset, tetra, fill=prev)) +
	geom_point(shape=21,alpha=0.5,size=3,colour='grey') + 
	scale_fill_gradient(high='red',low='blue') + 
	geom_hline(yintercept=0.5,linetype=2) +
	 theme_bw() +
        labs(x = 'Birth Strata', y='tetrachoric correlations', fill='Prevalence')
ggsave('Simulations_birth_year_dependent_prev_fix_assort.pdf')

save(df, file='Simulation_scene_1.RData')

# Scenario 2:
# Censoring occurred as a function of age as the issue of birth cohort
#    constant prevalence

assortvec <- c(0.1,0.2,0.3,0.4,0.5)
censorvec <- c(0.05,0.25,0.5) # scaling the censoring probability according to age
df <- data.frame()
for (assort in assortvec){
for (c in censorvec){
for (iiter in 1:100){
print(paste0(iiter))
agevec.1 <- runif(N, min = 20, max=90)
agevec.2 <- agevec.1 + runif(N, min = -5, max=5) # random assign the age gap, at most 10
prev.1 <- 0.2
prev.2 <- 0.2
liab.1 <- rnorm(N)
liab.2 <- (assort)*liab.1 + sqrt(1 - assort**2)*rnorm(N)
casevec.1 <- liab.1 > qnorm(prev.1, lower=F)
casevec.2 <- liab.2  > qnorm(prev.2, lower=F)
maskvec.1 <- (casevec.1 == T) & (runif(N) > c*(agevec.1/90))
maskvec.2 <- (casevec.2 == T) & (runif(N) > c*(agevec.2/90))
rvec0 <- c()
rvec1 <- c()
cen0 <- c()
for (ibin in 2:length(ageset)){
  inset <- (agevec.1 >= ageset[ibin-1]) & (agevec.1 < ageset[ibin])
  cen0 <- c(cen0, mean(c((casevec.1[inset] & !maskvec.1[inset]),(casevec.2[inset] & !maskvec.2[inset]))))
  rvec0 <- c(rvec0, cor(liab.1[inset],liab.2[inset]))
  rvec1 <- c(rvec1, tetrachoric(matrix(table(maskvec.1[inset],maskvec.2[inset]),2,2))$rho)
}
tmp <- data.frame(liab=rvec0, tetra=rvec1, cen=cen0, iter=iiter, assort=assort, age1 = ageset[1:(length(ageset)-1)], age2 = ageset[2:length(ageset)])
df <- rbind(df, tmp)
}}}

save(df, file='Simulation_scene_2.RData')

ggplot(df, aes(cen, tetra, fill=factor(assort))) +
        geom_point(shape=21,alpha=0.5,size=3,colour='grey') +
        geom_hline(yintercept=assortvec,linetype=2) +
        theme_bw() +
        labs(x = 'censoring', y='tetrachoric correlations', fill='Level of Assortment')
ggsave('Simulations_birth_year_dependent_censoring.pdf')

df$ageset <- paste0(df$age1, ' ~ ', df$age2)
ggplot(subset(df, assort==0.5), aes(ageset, tetra, fill=cen)) +
        geom_point(shape=21,alpha=0.5,size=3,colour='grey') +
        scale_fill_gradient(high='red',low='blue') +
        geom_hline(yintercept=0.5,linetype=2) +
         theme_bw() +
        labs(x = 'Birth Strata', y='tetrachoric correlations', fill='Censoring')
ggsave('Simulations_birth_year_dependent_censoring.fix_assort.pdf')













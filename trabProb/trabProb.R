f = function(x, b){
  (log(1+x)/(1+x)^(b+1))*((b^2*2^b)/(2^b-log(2)*b-1))
}
bGrid = seq(1, 10, 1)
col = rainbow(length(bGrid))
for (i in seq_along(bGrid)) {
  curve(f(x, bGrid[i]), add=i!=1, col=col[i],  ylim=c(0,4.3), 
        main = "Densidade para diferentes valores do parâmetro b", 
        xlab = "x", ylab = "f(x)")
}

legend(0.8,4.2,bGrid, col=col, lty=1, 
       title = "Valores para b", box.lty = 0)


Fx = function(x, b){
  -(2^b*(b*log(x+1)-(x+1)^b+1))/((2^b-log(2)*b-1)*(x+1)^b)
}

lik = function(b, x, n){
  lik = 2*log(b)*n+log(2)*b*n - log(2^b - log(2)*b - 1)*n + sum(log(log(x+1)))-b*sum(log(x+1))-sum(log(x+1))
  return(-lik)
}


for (i in seq_along(bGrid)) {
  curve(Fx(x, bGrid[i]), add=i!=1, col=col[i],  ylim=c(0,1),
        main = "Fda para diferentes valores do parâmetro b", 
        xlab = "x", ylab = "F(x)")
}

legend(0.875,0.8,bGrid, col=col, lty=1,
       title = "Valores para b", box.lty = 0)


library(GoFKernel)
bG = c(3,5,7,10)
ns = 1000
samp = matrix(0, ns, length(bG))
par = matrix(0, ns, length(bG))
vies = matrix(0, ns, length(bG))
eqm = matrix(0, ns, length(bG))

for (i in seq_along(bG)) {
  Fxx = function(x) Fx(x, b = bG[i])
  Finv = inverse(Fxx, lower = 0, upper = 1)
  for (k in 1:ns) {
    for (j in 1:ns) {
      samp[j,i] = Finv(runif(1))
    }
    par[k,i] = optim(0, lik, method = "L-BFGS-B", lower = 0.1, 
                    x = samp[,i], n = ns)$par
    vies[k,i] = par[k,i] - bG[i]
    eqm[k,i] = (par[k,i] - bG[i])^2
  }
}

df = data.frame("b" = bG, 
           "mean_bhat" = apply(par, 2, mean), 
           "var_bhat" = apply(par, 2, var),
           "vies" = apply(vies, 2, mean),
           "eqm" = apply(eqm, 2, mean))
xtable::xtable(df, digits = 3)

par(mfrow = c(2,2))

for (i in seq_along(bG)) {
  hist(par[,i], col = NULL, 
       border =  rainbow(length(bG))[i], 
       main = paste("Histograma do EMV b =", bG[i]), 
       xlab = expression(hat(b)), ylab = "Freq")
}


for (i in seq_along(bG)) {
  hist(samp[,i], probability = T, ylim = c(0,3.5), 
       col = NULL, border = "darkgray", 
       main = paste("Histograma para b =", bG[i]), 
       xlab = "x", ylab = "Densidade")
  lines(density(samp[,i]), add = T)
  curve(f(x, bG[i]), add=T, col= rainbow(length(bG))[i],
        ylim=c(0,3.5))
  legend(0.5,2.5, c("Densidade empírica", "Densidade teórica"), 
         col=c("black", rainbow(length(bG))[i]), lty=1,
         box.lty = 0)
}

library(readr)
library(readxl)
covid = read_delim("SESRSCOVID19_NOTIFICACOES.csv", 
                  ";", escape_double = FALSE, trim_ws = TRUE)
pop = read_excel("POP2020.xlsx")
covid$NM_MUNICIPIO = rm_accent(covid$NM_MUNICIPIO)
pop$`NOME DO MUNICÍPIO` = rm_accent(pop$`NOME DO MUNICÍPIO`)
covid$NM_MUNICIPIO = tolower(covid$NM_MUNICIPIO)
pop$`NOME DO MUNICÍPIO` = tolower(pop$`NOME DO MUNICÍPIO`)
names(covid)[2] = "muni"
names(pop)[4] = "muni"
df = merge(covid, pop, by = 'muni', all = T)
df$taxa = (df$POSITIVO_RTPCR+df$POSITIVO_TESTERAPIDO+df$POSITIVO_OUTROS)/df$`POPULAÇÃO ESTIMADA`

taxa = na.omit(df$taxa)
emv = optim(0, lik, method = "L-BFGS-B", lower = 0.000001, 
      x = taxa, n = length(taxa))$par

hist(taxa, probability = T, ylim = c(0,10), 
     col = NULL, border = "darkgray", 
     main = paste("Histograma incidência casos positivos covid RS 2020"), 
     xlab = "Proporção", ylab = "Densidade")
lines(density(taxa), add = T)
curve(f(x, emv), col= "purple",  ylim=c(0,8), add = T)
legend(0.2,8, c("Densidade empírica", "Densidade teórica EMV"), 
       col=c("black", "purple"), lty=1,
       box.lty = 0)

write.csv(df, "df.csv")

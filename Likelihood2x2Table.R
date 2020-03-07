table2x2.lik=function(table, lolim=0, hilim=100, main1 = 
     "Likelihood Functions For Odds Ratio\n",
     main2, OR2, dgt = 3, acc = 1, fignum)
{
###############################################################
## Function: 	orcon.lik  
## Author:		Richard Royall & Jeffrey Blume
## Origin:		November 1998 
## Revision:	February 2015
##
## Purpose:	Plot conditional likelihood function for odds 
##    		ratio in two-binomial model (condition on sum). 
##			Revision added profile likelihood, estimated 
##			likelihood, exact and pearson p-values, and
##			reports p-values based on maximum likelihood ratio 
##			as described in Choi et. al. 2015.
##
## Input:	x successes in m trials & y successes in n. 
##			table = c(a,b,c,d) 
##			where x=a, m=a+b and y=c, n=c+d 
##
## Notes:	Likelihood function evaluated at 1000*acc points
###############################################################

########################
#### Data/Preliminaries
########################

##	compute successes and trials in each group
	 x <- table[1,1]
	 m <- sum(table[1,])
	 y <- table[2,1]
	 n <- sum(table[2,])
 
########################
#### Hypothesis Testing
########################

## Fisher's exact Test
ft <- fisher.test(table)

## Pearson Chi-Square
ct <- chisq.test(table)

## Chi-square with continuity correction
ctc <- chisq.test(table,correct=FALSE)

## Get p-values 
pv <- data.frame(Fisher=ft$p.value,Pearson=ct$p.value,PearsonCC=ctc$p.value,row.names="p.value")

########################
#### Likelihoods
########################
 
## Conditional Likelihood (Non-Central Hypergeometric)

     k <- x + y
     s <- max(k - n, 0):min(m, k)
     lcoeff <- lgamma(x + 1) + lgamma(m - x + 1) + 
          lgamma(k - x + 1) + lgamma(n - k + x + 
          1) - lgamma(s + 1) - lgamma(m - s + 1) - 
          lgamma(k - s + 1) - lgamma(n - k + s + 
          1)
     if(x * y * (m - x) * (n - y) >= 1) {
          psi <- (x * (n - y))/(y * (m - x))
          std <- sqrt(1/x + 1/y + 1/(m - x) + 1/(
               n - y))
          ll <- exp(log(psi) - 5 * std)
          ul <- exp(log(psi) + 5 * std)
     }
     if(x * (n - y) == 0) {
          ll <- 0
          ul <- 2
     }
     if(y * (m - x) == 0) {
          ll <- 0
          ul <- 16
     }
     if(missing(lolim))
          lolim <- ll
     if(missing(hilim))
          hilim <- ul
     z <- c(1, if(missing(OR2)) 1 else OR2, seq(
          lolim, hilim, len = 1000 * acc))
     lik <- t(matrix(z, nrow = length(z), ncol = 
          length(s)))
     lik <- t(lik^(s - x))
     lik <- 1/(lik %*% exp(lcoeff))     

##   Now lik is the (unscaled) likelihood fn.

     if(y * (m - x) == 0) {
          lik <- lik * exp(lcoeff[length(s)])
          psihat <- NA
     }
     else {
          if(x * y * (m - x) * (n - y) >= 1 && {
               max(lik) == lik[length(lik)] || 
                    max(lik) == lik[3]
          }
          ) {
               llloc <- exp(log(psi) - 0.5 * 
                    std)
               ulloc <- exp(log(psi) + 0.5 * 
                    std)
               zloc <- seq(llloc, ulloc, len
                     = 500)
               likloc <- t(matrix(zloc, nrow
                     = length(zloc), ncol
                     = length(s)))
               likloc <- t(likloc^(s - x))
               likloc <- 1/(likloc %*% exp(
                    lcoeff))
               lik <- lik/max(likloc)
               psihat <- zloc[likloc == max(
                    likloc)]
          }
          else {
               lik <- lik/max(lik)
               psihat <- z[lik == max(lik)]
          }
     }
     if(missing(OR2))
          LR <- min(1/lik[1], 1000000)
     else LR <- min(lik[2]/lik[1], 1000000)
     z0 <- z[3:length(z)]
     lik0 <- lik[3:length(lik)]
     z1 <- z0[lik0 >= 1/8]
     i1 <- rep(1/8, length(z1))
     z2 <- z0[lik0 >= 1/32]
     i2 <- rep(1/32, length(z2))
     if(missing(main2)) {
          main2 <- paste("Observed proportions:", 
               x, "/", m, "and", y, "/", n)
          sub1 <- " "
     }
     else sub1 <- paste("Observed proportions:", x, 
               "/", m, "and", y, "/", n)
     
## Two-group Binomial Likelihood 
## The correct or 'true' likelihood is the produce of two binomials
## reparameterized from ( p1, p2 ) --> ( log(or), log(p2/(1-p2)) )
## This is just
## exp(psi*y1)*exp(lam*(y1+y2))*((1+exp(psi+lam))^(-n1))*((1+exp(lam))^(-n2))
##		where psi = log(OR) and lam = log(p2/(1-p2))

## Function to compute log-likelihood

loglik <- function(or,lam,tabf=table){log(or)*tabf[1,1]+
		lam*(sum(tabf[,1]))-sum(tabf[1,])*log(1+exp(log(or)+lam))-
		sum(tabf[2,])*log(1+exp(lam))}

## Set range of oddsratio and logOdds grp2 (hidden)
lo.or = 0 ; hi.or = hilim
lo.lam = -15 ; hi.lam = 15

or = seq(lo.or,hi.or,len=(10000*acc+1))
lam = seq(lo.lam,hi.lam,len=(1000*acc+1))

## Find the global MLE for lambda based on group 2 ; reparameterized [pg 374]
lam.max = log(table[2,1]/table[2,2])	
                                    
or.null = which.min(abs(or - 1))   	
	# instead of which(or==1) b/c 1 might not be in the or vector

## Estimated likelihood (uses lam.max for plug-in of nusiance parameter lambda)
	# Note: The estimted likelihood essentially assumes lam.max is 'correct'  
	#		so the likelihood is too concentrated because estimate is
	#		treated as a known constant.
llik.est = loglik(or,lam.max)		         
lik.est  = exp(llik.est-max(llik.est))    # scale the function between 0 and 1

lre = 1/lik.est[or.null]       
	# The numerator is max(lik.est)==1 because of the scaling 
	# this is the LR for the MLE versus 1 (for the or)

pv.e = 1-pchisq((2*log(lre)),df=1)        
	# -2*log(LR) ~ chisq(df=1) ; forumla from Choi et. al. (2015) [pg 17] 

lr.e = data.frame(p.value=pv.e,maxLR=lre,row.names="Est.Like") 

## Profile likelihood (uses the MLE lambda for each OR )
	# Note: The profile likelihood plugs in the best estimator for lambda 
	#		conditional on each OR. That is, for each OR, it finds the best 
	#		estimate of lambda. As a result, the plug-in estimate of lambda 
	#		changes depending on the OR. 
llik.pro = outer(or,lam,"loglik",tabf=table) 	# "order" tries all combinations
llik.pro = apply(llik.pro,1,max)		# maximizes likelihood over lambda for each OR
lik.pro  = exp(llik.pro-max(llik.pro)) 	# scale the function between 0 and 1

lrp = 1/lik.pro[or.null] 
	# The numerator is max(lik.pro)==1 because of the scaling 
	# this is the LR for the MLE versus 1 (for the or)

pv.p = 1-pchisq((2*log(lrp)),df=1)
	# -2*log(LR) ~ chisq(df=1) ; forumla from Choi et. al. (2015) [pg 17] 

lr.p = data.frame(p.value=pv.p,maxLR=lrp,row.names="Prof.Like") 

########################
#### Plotting
########################

	 plot(z0, lik0, xlim = c(lolim, hilim), ylim = 
          range(lik0, 0, 1), type = "l", main = 
          paste(if(!missing(fignum)) fignum, 
          main1, main2), xlab = "Odds Ratio", 
          ylab = " ", sub = sub1)
     lines(z1, i1, type = "l")
     lines(z2, i2, type = "l")
	 
	 lines(or,lik.est,col='forestgreen') 	# estimated likelihood
	 lines(or,lik.pro,col='royalblue')		# Profile likelihood
	 abline(v=1,lty=2,lwd=0.5)				# Null line
	 
	 legend(3*(hi.or-lo.or)/4,0.8,c("Conditional","Estimated","Profile"),
	 		lty=1,col=c("black","forestgreen","royalblue"),bty="n")
	 
     whr <- if(psihat == "NA") (4 * lolim + hilim)/5
           else if(psihat < (lolim + hilim)/2)
          (lolim + 4 * hilim)/5
     else (4 * lolim + hilim)/5
     if(!missing(fignum))
          text(whr, 0.87, paste("L(", OR2, 
               ")/L(1)=", signif(LR, dgt)))
     else {
          text(whr, 0.95, paste("Max at", signif(
               psihat, 3), "\n1/8 LI (", if(
               min(z1) == z[3] && z[3] == 0) 
                    "0" else if(min(z1) == 
               z[3])
               "NA"
          else if(min(z1) == "NA")
               "NA"
          else round(min(z1), 2), ",", if(max(z1) ==
               z[length(z)]) "NA" else if(max(
               z1) == "NA")
               "NA"
          else round(max(z1), 2), ")", 
               "\n1/32 LI (", if(min(z2) == z[
               3] && z[3] == 0) "0" else if(
               min(z2) == z[3])
               "NA"
          else if(min(z2) == "NA")
               "NA"
          else round(min(z2), 2), ",", if(max(z2) ==
               z[length(z)]) "NA" else if(max(
               z2) == "NA")
               "NA"
          else round(max(z2), 2), ")", "\nL(", 
               if(missing(OR2)) signif(psihat, 
                    dgt) else signif(OR2, 
                    dgt), if(LR == 1000000
               ) ")/L(1)>" else ")/L(1)=", 
               signif(LR, dgt)))
     }

## Show p-values and conditional LR p-value
## as suggested in Choi et. al. (2015)
pv.c = 1-pchisq((2*log(LR)),df=1)
	# -2*log(LR) ~ chisq(df=1) ; forumla from Choi et. al. (2015) [pg 17] 

lr.c = data.frame(p.value=pv.c,maxLR=LR,row.names="Cond.Like") 

out=rbind(data.frame(cbind(t(pv),maxLR=rep(NA,3))),lr.c,lr.e,lr.p)
print(out, digits=5)	 
	 
}

#### Example 1 ####
table1=matrix(c(1,9,5,5),ncol=2,byrow=TRUE)   # Example from Choi et. al. (2015)
table2x2.lik(table1,hilim=2)

#### Example 2 ####
table2=matrix(c(11,3,1,9),ncol=2,byrow=TRUE)
table2x2.lik(table2,hilim=100)

#### Example 1 ####
table3=matrix(c(16,14,12,68),ncol=2,byrow=TRUE)
table2x2.lik(table3,hilim=30)


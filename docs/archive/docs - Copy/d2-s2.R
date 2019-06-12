
library(EpiModel)
mynet <- network.initialize(256, directed=FALSE)
age <- rep(18:25, each=32)
mynet %v%'age' <- age

formation <- ~edges+concurrent+absdiff("age")+ degrange(from = 4)
target.stats <- c(152, 80, 192, 0)

myfit <- netest(mynet,
            formation=formation,
            target.stats = target.stats,
            coef.diss = dissolution_coefs(~offset(edges), 90))

mydx <- netdx(myfit, nsims=10, nsteps=100)
mydx
boxplot(mydx$stats[[1]])

mycontrol <- control.net("SIS", nsteps = 100, nsims = 10,
                verbose = TRUE)
myinit <- init.net(i.num = 10)
myparam <-param.net(inf.prob = 0.5, act.rate = 0.6,
                rec.rate = 0.1)
mySIS <- netsim(myfit, param = myparam, control = mycontrol,
                init = myinit)
plot(mySIS)

mySIS$stats$nwstats
plot(mySIS, type = "formation", sim.lines = TRUE)

race <- rep(c("B","W"), 128)
xtabs(~race+age)
mynet %v%'race' <- race
formation <- ~edges+concurrent+absdiff("age")+degrange(from = 4)+
                nodefactor("race")+nodematch("race")
target.stats <- c(152, 80, 192, 0, 176, 120)
myfit <- netest(mynet,
                formation=formation,
                target.stats = target.stats,
                coef.diss = dissolution_coefs(~offset(edges), 90))

mydx <- netdx(myfit, nsims=10, nsteps=100)
mydx
boxplot(mydx$stats[[1]])


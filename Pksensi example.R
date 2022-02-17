library(pksensi)
library(httk)


pbtk1cpt <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dAgutlument = - kgutabs * Agutlument
    dAcompartment = kgutabs * Agutlument - ke * Acompartment
    dAmetabolized = ke * Acompartment
    Ccompartment = Acompartment / vdist * BW;
    list(c(dAgutlument, dAcompartment, dAmetabolized), 
         "Ccompartment" = Ccompartment) 
  })
}

initState <- c(Agutlument = 10, Acompartment = 0, Ametabolized = 0)
t <- seq(from = 0.01, to = 24.01, by = 1)
outputs <- c("Ccompartment")


pars1comp <- (parameterize_1comp(chem.name = "acetaminophen"))
#> Human volume of distribution returned in units of L/kg BW.

q <- c("qunif", "qunif", "qunif", "qnorm")
q.arg <- list(list(min = pars1comp$Vdist / 2, max = pars1comp$Vdist * 2),
              list(min = pars1comp$kelim / 2, max = pars1comp$kelim * 2),
              list(min = pars1comp$kgutabs / 2, max = pars1comp$kgutabs * 2),
              list(mean = pars1comp$BW, sd = 5))

set.seed(1234)
params <- c("vdist", "ke", "kgutabs", "BW")
x <- rfast99(params, n = 200, q = q, q.arg = q.arg, replicate = 1)

out <- solve_fun(x, time = t, func = pbtk1cpt, initState = initState, outnames = outputs)

pksim(out) 

plot(out)

@defcomp biodiversity begin
    # Number of species
    nospecies = Variable(index=[time])

    # additive parameter
    bioloss = Parameter(default=0.003)

    # multiplicative parameter
    biosens = Parameter(default=0.001)

    # Temperature
    temp = Parameter(index=[time])

    # benchmark temperature change
    dbsta = Parameter()

    # Number of species in the year 2000
    nospecbase = Parameter()

    function run_timestep(p, v, d, t)
        
        if gettime(t) > 2000
            dt = abs(p.temp[t] - p.temp[t - 1])

            v.nospecies[t] = max(
            p.nospecbase / 100,
            v.nospecies[t - 1] * (1.0 - p.bioloss - p.biosens * dt * dt / p.dbsta / p.dbsta)
            )
        else
            v.nospecies[t] = p.nospecbase
        end
    end
end

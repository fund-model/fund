@defcomp impacttropicalstorms begin
    regions = Index()

    hurrdam = Variable(index=[time,regions])
    hurrdead = Variable(index=[time,regions])

    hurrbasedam = Parameter(index=[regions])
    hurrdamel = Parameter(default=-0.514)
    hurrnl = Parameter(default=3)
    hurrpar = Parameter(default=0.04)
    hurrbasedead = Parameter(index=[regions])
    hurrdeadel = Parameter(default=-0.501)

    gdp90 = Parameter(index=[regions])
    pop90 = Parameter(index=[regions])

    population = Parameter(index=[time,regions])
    income = Parameter(index=[time,regions])

    regstmp = Parameter(index=[time,regions])

    function run_timestep(p, v, d, t)

        for r in d.regions
            ypc = p.income[t, r] / p.population[t, r] * 1000.0
            ypc90 = p.gdp90[r] / p.pop90[r] * 1000.0

            # This is hurrican damage
            v.hurrdam[t, r] = 0.001 * p.hurrbasedam[r] * p.income[t, r] * (ypc / ypc90)^p.hurrdamel * ((1.0 + p.hurrpar * p.regstmp[t, r])^p.hurrnl - 1.0)

            v.hurrdead[t, r] = 1000.0 * p.hurrbasedead[r] * p.population[t, r] * (ypc / ypc90)^p.hurrdeadel * ((1.0 + p.hurrpar * p.regstmp[t, r])^p.hurrnl - 1.0)
        end
    end
end

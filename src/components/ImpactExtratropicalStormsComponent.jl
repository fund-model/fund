@defcomp impactextratropicalstorms begin
    regions = Index()

    extratropicalstormsdam = Variable(index=[time,regions])
    extratropicalstormsdead = Variable(index=[time,regions])

    extratropicalstormsbasedam = Parameter(index=[regions])
    extratropicalstormsdamel = Parameter(default=-0.514)
    extratropicalstormspar = Parameter(index=[regions])
    extratropicalstormsbasedead = Parameter(index=[regions])
    extratropicalstormsdeadel = Parameter(default=-0.501)
    extratropicalstormsnl = Parameter(default=1)

    gdp90 = Parameter(index=[regions])
    pop90 = Parameter(index=[regions])

    population = Parameter(index=[time,regions])
    income = Parameter(index=[time,regions])

    acco2 = Parameter(index=[time])
    co2pre = Parameter()

    function run_timestep(p, v, d, t)

        for r in d.regions
            ypc = p.income[t, r] / p.population[t, r] * 1000.0
            ypc90 = p.gdp90[r] / p.pop90[r] * 1000.0

            v.extratropicalstormsdam[t, r] = p.extratropicalstormsbasedam[r] * p.income[t, r] * (ypc / ypc90)^p.extratropicalstormsdamel * ((1.0 + (p.extratropicalstormspar[r] * (p.acco2[t] / p.co2pre)))^p.extratropicalstormsnl - 1.0)
            v.extratropicalstormsdead[t, r] = 1000.0 * p.extratropicalstormsbasedead[r] * p.population[t, r] * (ypc / ypc90)^p.extratropicalstormsdeadel * ((1.0 + (p.extratropicalstormspar[r] * (p.acco2[t] / p.co2pre)))^p.extratropicalstormsnl - 1.0)
        end
    end
end

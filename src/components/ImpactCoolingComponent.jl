@defcomp impactcooling begin
    regions = Index()

    cooling = Variable(index=[time,regions])

    cebm = Parameter(index=[regions])
    gdp90 = Parameter(index=[regions])

    population = Parameter(index=[time,regions])
    pop90 = Parameter(index=[regions])

    income = Parameter(index=[time,regions])
    ceel = Parameter(default=0.8)

    temp = Parameter(index=[time,regions])
    cenl = Parameter(default=1.5)

    cumaeei = Parameter(index=[time,regions])

    function run_timestep(p, v, d, t)

        if !is_first(t)
            for r in d.regions
                ypc = p.income[t, r] / p.population[t, r] * 1000.0
                ypc90 = p.gdp90[r] / p.pop90[r] * 1000.0

                v.cooling[t, r] = p.cebm[r] * p.cumaeei[t, r] * p.gdp90[r] * (p.temp[t, r] / 1.0)^p.cenl * (ypc / ypc90)^p.ceel * p.population[t, r] / p.pop90[r]
            end
        end
    end
end

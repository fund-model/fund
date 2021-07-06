@defcomp impactheating begin
    regions = Index()

    heating = Variable(index=[time,regions])

    hebm = Parameter(index=[regions])
    gdp90 = Parameter(index=[regions])

    population = Parameter(index=[time,regions])
    pop90 = Parameter(index=[regions])

    income = Parameter(index=[time,regions])
    heel = Parameter(default=0.8)

    temp = Parameter(index=[time,regions])
    henl = Parameter(default=0.5)

    cumaeei = Parameter(index=[time,regions])

    function run_timestep(p, v, d, t)

        if !is_first(t)
            for r in d.regions
                ypc = p.income[t, r] / p.population[t, r] * 1000.0
                ypc90 = p.gdp90[r] / p.pop90[r] * 1000.0

                v.heating[t, r] = p.hebm[r] * p.cumaeei[t, r] * p.gdp90[r] * atan(p.temp[t, r]) / atan(1.0) * (ypc / ypc90)^p.heel * p.population[t, r] / p.pop90[r]
            end
        end
    end
end
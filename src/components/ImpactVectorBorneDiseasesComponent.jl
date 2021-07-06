@defcomp impactvectorbornediseases begin
    regions = Index()

    dengue = Variable(index=[time,regions])
    schisto = Variable(index=[time,regions])
    malaria = Variable(index=[time,regions])

    dfbs = Parameter(index=[regions])
    dfch = Parameter(index=[regions])
    smbs = Parameter(index=[regions])
    smch = Parameter(index=[regions])
    malbs = Parameter(index=[regions])
    malch = Parameter(index=[regions])

    gdp90 = Parameter(index=[regions])
    pop90 = Parameter(index=[regions])

    income = Parameter(index=[time,regions])
    population = Parameter(index=[time,regions])
    temp = Parameter(index=[time,regions])

    dfnl = Parameter(default=1)
    vbel = Parameter(default=-2.65)
    smnl = Parameter(default=1)
    malnl = Parameter(default=1)

    function run_timestep(p, v, d, t)

        for r in d.regions
            ypc = 1000.0 * p.income[t, r] / p.population[t, r]
            ypc90 = p.gdp90[r] / p.pop90[r] * 1000.0

            v.dengue[t, r] = p.dfbs[r] * p.population[t, r] * p.dfch[r] * p.temp[t, r]^p.dfnl * (ypc / ypc90)^p.vbel

            v.schisto[t, r] = p.smbs[r] * p.population[t, r] * p.smch[r] * p.temp[t, r]^p.smnl * (ypc / ypc90)^p.vbel

            if v.schisto[t, r] < -p.smbs[r] * p.population[t, r] * (ypc / ypc90)^p.vbel
                v.schisto[t, r] = -p.smbs[r] * p.population[t, r] * (ypc / ypc90)^p.vbel
            end

            v.malaria[t, r] = p.malbs[r] * p.population[t, r] * p.malch[r] * p.temp[t, r]^p.malnl * (ypc / ypc90)^p.vbel
        end
    end
end
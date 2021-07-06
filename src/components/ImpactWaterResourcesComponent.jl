@defcomp impactwaterresources begin
    regions = Index()

    watech = Variable(index=[time])
    watechrate = Parameter(default=0.005)

    water = Variable(index=[time,regions])
    wrbm = Parameter(index=[regions])
    wrel = Parameter(default=0.85)
    wrnl = Parameter(default=1)
    wrpl = Parameter(default=0.85)

    gdp90 = Parameter(index=[regions])
    income = Parameter(index=[time,regions])

    population = Parameter(index=[time,regions])
    pop90 = Parameter(index=[regions])

    temp = Parameter(index=[time,regions])

    function run_timestep(p, v, d, t)

        if gettime(t) > 2000
            v.watech[t] = (1.0 - p.watechrate)^(gettime(t) - 2000)
        else
            v.watech[t] = 1.0
        end

        for r in d.regions
            ypc = p.income[t, r] / p.population[t, r] * 1000.0
            ypc90 = p.gdp90[r] / p.pop90[r] * 1000.0

            water = p.wrbm[r] * p.gdp90[r] * v.watech[t] * (ypc / ypc90)^p.wrel * (p.population[t, r] / p.pop90[r])^p.wrpl * p.temp[t, r]^p.wrnl

            if water > 0.1 * p.income[t, r]
                v.water[t, r] = 0.1 * p.income[t, r]
            else
                v.water[t, r] = water
            end
        end
    end
end
@defcomp impactagriculture begin
    regions = Index()

    gdp90 = Parameter(index=[regions])
    income = Parameter(index=[time,regions])
    pop90 = Parameter(index=[regions])
    population = Parameter(index=[time,regions])

    agrish = Variable(index=[time,regions])
    agrish0 = Parameter(index=[regions])
    agel = Parameter(default=0.31)

    agrate = Variable(index=[time,regions])
    aglevel = Variable(index=[time,regions])
    agco2 = Variable(index=[time,regions])
    agcost = Variable(index=[time,regions])

    agrbm = Parameter(index=[regions])
    agtime = Parameter(index=[regions])
    agnl = Parameter(default=2)

    aglparl = Parameter(index=[regions])
    aglparq = Parameter(index=[regions])

    agcbm = Parameter(index=[regions])
    co2pre = Parameter()

    temp = Parameter(index=[time,regions])
    acco2 = Parameter(index=[time])

    function run_timestep(p, v, d, t)
        
        DBsT = 0.04     # base case yearly warming
        
        if is_first(t)
            for r in d.regions
                v.agrate[t, r] = p.agrbm[r] * (0.005 / DBsT)^p.agnl * p.agtime[r]
            end
        else
            for r in d.regions
                ypc = p.income[t, r] / p.population[t, r] * 1000.0
                ypc90 = p.gdp90[r] / p.pop90[r] * 1000.0

                v.agrish[t, r] = p.agrish0[r] * (ypc / ypc90)^(-p.agel)
            end

            for r in d.regions
                dtemp = abs(p.temp[t, r] - p.temp[t - 1, r])

                if isnan((dtemp / 0.04)^p.agnl)
                    v.agrate[t, r] = 0.0
                else
                    v.agrate[t, r] = p.agrbm[r] * (dtemp / 0.04)^p.agnl + (1.0 - 1.0 / p.agtime[r]) * v.agrate[t - 1, r]
                end
            end

            for r in d.regions
                v.aglevel[t, r] = p.aglparl[r] * p.temp[t, r] + p.aglparq[r] * p.temp[t, r]^2.0
            end

            for r in d.regions
                v.agco2[t, r] = p.agcbm[r] / log(2.0) * log(p.acco2[t - 1] / p.co2pre)
            end

            for r in d.regions
                v.agcost[t, r] = min(1.0, v.agrate[t, r] + v.aglevel[t, r] + v.agco2[t, r]) * v.agrish[t, r] * p.income[t, r]
            end
        end
    end
end

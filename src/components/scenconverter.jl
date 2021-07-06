@defcomp scenconverter begin
    regions = Index()

    population = Parameter(index=[time,regions])
    income = Parameter(index=[time,regions])
    energuse = Parameter(index=[time,regions])
    emission = Parameter(index=[time,regions])

    pop0 = Variable(index=[regions])
    scenpgrowth = Variable(index=[time,regions])
    gdp0 = Variable(index=[regions])
    scenypcgrowth = Variable(index=[time,regions])
    scenaeei = Variable(index=[time,regions])
    scenacei = Variable(index=[time,regions])
    energint0 = Variable(index=[regions])
    emissint0 = Variable(index=[regions])
    energint = Variable(index=[time,regions])
    emissint = Variable(index=[time,regions])

    function run_timestep(p, v, d, t)

        if is_first(t)
            for r in d.regions
                v.pop0[r] = p.population[t, r]
                v.gdp0[r] = p.income[t, r]
                v.energint0[r] = p.energuse[t,r] / p.income[t, r]
                v.emissint0[r]  = p.emission[t, r] / p.energuse[t, r]

                v.energint[t, r] = v.energint0[r]
                v.emissint[t, r] = v.emissint0[r]
            end
        end

        for r in d.regions
            if t.t < 1050
                v.scenpgrowth[t, r] = (p.population[t + 1, r] / p.population[t, r] - 1.) * 100.
                v.scenypcgrowth[t, r] = (p.income[t + 1, r] / p.income[t, r] / (1 + 0.01 * v.scenpgrowth[t, r]) - 1.) * 100.
            end

            v.energint[t, r] = p.energuse[t, r] / p.income[t,r]
            v.emissint[t, r] = p.emission[t, r] / p.energuse[t, r]

            if !is_first(t)
                v.scenaeei[t, r] = -(v.energint[t, r] / v.energint[t - 1, r] - 1.) * 100.
                v.scenacei[t, r] = -(v.emissint[t, r] / v.emissint[t - 1, r] - 1.) * 100.
            end
        end
    end
end

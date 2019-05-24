func_U(c, γ) = γ==1.0 ? log(c) : c^(1-γ)/(1-γ)
func_inv_U(u, γ) = γ==1.0 ? exp(u) : (u*(1-γ))^(1/(1-γ))
func_V(c, η) = η==1.0 ? log(c) : c^(1-η)/(1-η)

@defcomp welfare begin
    regions = Index()

    prtp = Parameter(default=0.015)
    temporal_inequality_aversion = Parameter(default=1.5)
    regional_inequality_aversion = Parameter(default=0.)

    start_timeperiod::Int = Parameter(default=2020)

    populationin1 = Parameter(index=[time,regions])
    consumption = Parameter(index=[time,regions])

    c_ede = Variable(index=[time])

    welfare = Variable(index=[time])
    marginal_welfare = Variable(index=[time, regions])
    cum_welfare = Variable(index=[time])


    function run_timestep(p, v, d, t)
        year = gettime(t)

        if year >= p.start_timeperiod
            global_population = sum(p.populationin1[t,r] for r in d.regions)

            c_ede = 0.
            for r in d.regions
                c_ede += func_U(p.consumption[t,r] / p.populationin1[t,r], p.regional_inequality_aversion) * p.populationin1[t,r] / global_population
            end
            v.c_ede[t] = func_inv_U(c_ede, p.regional_inequality_aversion)

            v.welfare[t] = global_population * func_V(v.c_ede[t], p.temporal_inequality_aversion) * 1.0/(1+p.prtp)^(year-p.start_timeperiod)

            if year==p.start_timeperiod
                v.cum_welfare[t] = v.welfare[t]
            else
                v.cum_welfare[t] = v.cum_welfare[t-1] + v.welfare[t]
            end

            for r in d.regions
                v.marginal_welfare[t,r] = v.c_ede[t]^(p.regional_inequality_aversion-p.temporal_inequality_aversion) * (p.consumption[t,r] / p.populationin1[t,r])^(-p.regional_inequality_aversion) * 1.0/(1+p.prtp)^(year-p.start_timeperiod)
            end            
        end
    end
end

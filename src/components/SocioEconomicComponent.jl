@defcomp socioeconomic begin
	regions = Index()

	income = Variable(index=[time,regions])
	consumption = Variable(index=[time,regions])
	ypc = Variable(index=[time,regions])
	ygrowth = Variable(index=[time,regions])

	plus = Variable(index=[time,regions])
	urbpop = Variable(index=[time,regions])
	popdens = Variable(index=[time,regions])

	globalconsumption = Variable(index=[time])
	globalypc = Variable(index=[time])
	globalincome = Variable(index=[time])

	ypc90 = Variable(index=[regions])


	pgrowth             = Parameter(index=[time,regions])
	ypcgrowth           = Parameter(index=[time,regions])
	eloss               = Parameter(index=[time,regions])
	sloss               = Parameter(index=[time,regions])
	mitigationcost      = Parameter(index=[time,regions])
	area                = Parameter(index=[time,regions])
	globalpopulation    = Parameter(index=[time])
	population          = Parameter(index=[time,regions])
	populationin1       = Parameter(index=[time,regions])

	plus90  = Parameter(index=[regions])
	gdp90   = Parameter(index=[regions])
	pop90   = Parameter(index=[regions])
	urbcorr = Parameter(index=[regions])
	gdp0    = Parameter(index=[regions])


	runwithoutdamage = Parameter{Bool}(default=false)
	consleak    = Parameter(default=0.25)
	plusel      = Parameter(default=0.25)
    savingsrate = Parameter(default=0.2)

    function run_timestep(p, v, d, t)

        if is_first(t)
            for r in d.regions
                v.income[t, r] = p.gdp0[r]
                v.ypc[t, r] = v.income[t, r] / p.population[t, r] * 1000.0
                v.consumption[t, r] = v.income[t, r] * 1000000000.0 * (1.0 - p.savingsrate)
            end

            v.globalconsumption[t] = sum(v.consumption[t,:])

            for r in d.regions
                v.ypc90[r] = p.gdp90[r] / p.pop90[r] * 1000
            end
        else

            # Calculate income growth rate
            for r in d.regions
                v.ygrowth[t, r] = (1 + 0.01 * p.pgrowth[t - 1, r]) * (1 + 0.01 * p.ypcgrowth[t - 1, r]) - 1
            end

            # Calculate income
            for r in d.regions
                oldincome = v.income[t - 1, r] - (gettime(t) >= 1990 && !p.runwithoutdamage ? p.consleak * p.eloss[t - 1, r] / 10.0 : 0)

                v.income[t, r] = (1 + v.ygrowth[t, r]) * oldincome - p.mitigationcost[t - 1, r]
            end

            # Check for unrealistic values
            for r in d.regions
                if v.income[t, r] < 0.01 * p.population[t, r]
                    v.income[t, r] = 0.1 * p.population[t, r]
                end
            end

            for r in d.regions
                v.ypc[t, r] = v.income[t, r] / p.population[t, r] * 1000.0
            end

            for r in d.regions
                v.consumption[t, r] = max(v.income[t, r] * 1000000000.0 * (1.0 - p.savingsrate) - (p.runwithoutdamage ? 0.0 :   (p.eloss[t - 1, r] + p.sloss[t - 1, r]) * 1000000000.0), 0.0)
            end
            v.globalconsumption[t] = sum(v.consumption[t,:])

            for r in d.regions
                v.plus[t, r] = p.plus90[r] * (v.ypc[t, r] / v.ypc90[r])^p.plusel

                if v.plus[t, r] > 1
                    v.plus[t, r] = 1.0
                end
            end

            for r in d.regions
                v.popdens[t, r] = p.population[t, r] / p.area[t, r] * 1000000.0
            end

            for r in d.regions
                # ERROR This doesn't make sense for t < 40
                v.urbpop[t, r] = (0.031 * sqrt(v.ypc[t, r]) - 0.011 * sqrt(v.popdens[t, r])) / (1.0 + 0.031 * sqrt(v.ypc[t,r]) - 0.011 * sqrt(v.popdens[t, r])) / (1 + p.urbcorr[r] / (1 + 0.001 * (gettime(t) - 1990)^2.))
                # DA: urbcorr needs to be changed to a function if this is to be made uncertain
            end
            
            v.globalincome[t] = sum(v.income[t,:])

            v.globalypc[t] = sum(v.income[t,:] .* 1000000000.0) / sum(p.populationin1[t,:])
        end
    end
end
@defcomp impactsealevelrise begin
    regions = Index()

    wetval = Variable(index=[time,regions]) # Valuation of wetlands
    wetlandloss = Variable(index=[time,regions]) # Amount of wetland lost
    cumwetlandloss = Variable(index=[time,regions]) # Cumulative wetland loss over time
    wetlandgrowth = Variable(index=[time,regions])
    wetcost = Variable(index=[time,regions]) # The cost of wetland actually lost in this year


    dryval = Variable(index=[time,regions]) # Valuation of dryland
    landloss = Variable(index=[time,regions]) # Dryland lost in this year
    cumlandloss = Variable(index=[time,regions]) # Cumulative dryland loss over time
    drycost = Variable(index=[time,regions]) # The cost of dryland actually lost in this year

    npprotcost = Variable(index=[time,regions])
    npwetcost = Variable(index=[time,regions])
    npdrycost = Variable(index=[time,regions])

    protlev = Variable(index=[time,regions]) # Share of coastline protected in this time period
    protcost = Variable(index=[time,regions]) # The cost of protection in this time period

    enter = Variable(index=[time,regions]) # Migration in (in people)
    leave = Variable(index=[time,regions]) # Migration out (in people)
    entercost = Variable(index=[time,regions]) # Cost to the economy of immigration
    leavecost = Variable(index=[time,regions]) # Cost to the economy of emigration

    imigrate = Variable(index=[regions,regions]) # Migration matrix

    incdens     = Parameter(default=0.000635) # Normalization income density
    emcst       = Parameter(default=3) # emigration loss benchmark value
    immcst      = Parameter(default=0.4) # immigration loss benchmark
    dvydl       = Parameter(default=1) # income density elasticity of dryland value
    wvel        = Parameter(default=1.16) # wetland value income elasticity
    wvbm        = Parameter(default=0.00588) # wetland value benchmark value
    slrwvpopdens0 = Parameter(default=27.5937717888728) # wetland value income normalization value
    wvpdl       = Parameter(default=0.47)
    wvsl        = Parameter(default=-0.11)
    dvbm        = Parameter(default=0.004) # dryland value benchmark value
    slrwvypc0   = Parameter(default=25000)

    pc          = Parameter(index=[regions])
    slrprtp     = Parameter(index=[regions])
    wmbm        = Parameter(index=[regions])
    dlbm        = Parameter(index=[regions]) # dryland loss benchmark value
    drylandlossparam = Parameter(index=[regions])
    wlbm        = Parameter(index=[regions]) # wetland loss benchmark per slr change
    coastpd     = Parameter(index=[regions]) # coastal population density
    wetmax      = Parameter(index=[regions])
    wetland90   = Parameter(index=[regions])
    maxlandloss = Parameter(index=[regions])

    sea = Parameter(index=[time]) # sea-level rise in m since pre-industrial

    migrate = Parameter(index=[regions,regions])

    income = Parameter(index=[time,regions])
    population = Parameter(index=[time,regions])
    area = Parameter(index=[time,regions]) # Land area km^2

    function run_timestep(p, v, d, t)

        if is_first(t)
            for r1 in d.regions
                for r2 in d.regions
                    immsumm = 0
                    for i in d.regions
                        immsumm += p.migrate[i, r1]
                    end
                    v.imigrate[r1, r2] = p.migrate[r2, r1] / immsumm
                end

                t0 = TimestepIndex(1)
                v.landloss[t0, r1] = 0.0
                v.cumlandloss[t0, r1] = 0.0
                v.cumwetlandloss[t0, r1] = 0.0
                v.wetlandgrowth[t0, r1] = 0.0
            end
        else
            # slr in m/year
            ds = p.sea[t] - p.sea[t - 1]

            for r in d.regions
                ypc = p.income[t, r] / p.population[t, r] * 1000.0 # income per capita
                ypcprev = p.income[t - 1, r] / p.population[t - 1, r] * 1000.0
                ypcgrowth = ypc / ypcprev - 1.0

                if gettime(t) == 1951
                    ypcgrowth = 0
                end

                # Needs to be in $bn per km^2
                # Income is in billion, area is in km^2
                incomedens = p.income[t, r] / p.area[t, r]

                incomedensprev = p.income[t - 1, r] / p.area[t - 1, r]

                incomedensgrowth = incomedens / incomedensprev - 1.0

                # In population/km^2
                ## population is in million, area is in km^2
                popdens = p.population[t, r] / p.area[t, r] * 1000000.0
                popdensprev = p.population[t - 1, r] / p.area[t - 1, r] * 1000000.0
                popdensgrowth = popdens / popdensprev - 1.0

                # Unit of dryval is $bn/km^2
                v.dryval[t, r] = p.dvbm * (incomedens / p.incdens)^p.dvydl

                # Unit of wetval is $bn/km^2
                v.wetval[t, r] = p.wvbm *
                    (ypc / p.slrwvypc0)^p.wvel *
                    (popdens / p.slrwvpopdens0)^p.wvpdl *
                    ((p.wetland90[r] - v.cumwetlandloss[t - 1, r]) / p.wetland90[r])^p.wvsl

                potCumLandloss = min(p.maxlandloss[r], p.dlbm[r] * p.sea[t]^p.drylandlossparam[r])

                potLandloss = potCumLandloss - v.cumlandloss[t - 1, r]

                ## If sea levels fall, no protection is build
                if ds < 0
                    v.npprotcost[t, r] = 0
                    v.npwetcost[t, r] = 0
                    v.npdrycost[t, r] = 0
                    v.protlev[t, r] = 0
                # If the discount rate < -100% people will not build protection
                elseif (1.0 + p.slrprtp[r] + ypcgrowth) < 0.0
                    v.npprotcost[t, r] = 0
                    v.npwetcost[t, r] = 0
                    v.npdrycost[t, r] = 0
                    v.protlev[t, r] = 0
                # Dryland value is worthless
                elseif (1.0 + p.dvydl * incomedensgrowth) < 0.0
                    v.npprotcost[t, r] = 0
                    v.npwetcost[t, r] = 0
                    v.npdrycost[t, r] = 0
                    v.protlev[t, r] = 0
                # Is protecting the coast infinitly expensive?
                elseif (1.0 / (1.0 + p.slrprtp[r] + ypcgrowth)) >= 1
                    v.npprotcost[t, r] = 0
                    v.npwetcost[t, r] = 0
                    v.npdrycost[t, r] = 0
                    v.protlev[t, r] = 0
                # Is dryland infinitly valuable?
                elseif ((1.0 + p.dvydl * incomedensgrowth) / (1.0 + p.slrprtp[r] + ypcgrowth)) >= 1.0
                    v.npprotcost[t, r] = 0
                    v.npwetcost[t, r] = 0
                    v.npdrycost[t, r] = 0
                    v.protlev[t, r] = 1
                # Is wetland infinitly valuable?
                elseif ((1.0 + p.wvel * ypcgrowth + p.wvpdl * popdensgrowth + p.wvsl * v.wetlandgrowth[t - 1, r]) / (1.0 + p.slrprtp[r] + ypcgrowth)) >= 1.0
                    v.npprotcost[t, r] = 0
                    v.npwetcost[t, r] = 0
                    v.npdrycost[t, r] = 0
                    v.protlev[t, r] = 0
                else
                    # NPV of protecting the whole coast
                    # pc is in $bn/m
                    v.npprotcost[t, r] = p.pc[r] * ds * (1.0 + p.slrprtp[r] + ypcgrowth) / (p.slrprtp[r] + ypcgrowth)

                    # NPV of wetland
                    if (1.0 + p.wvel * ypcgrowth + p.wvpdl * popdensgrowth + p.wvsl * v.wetlandgrowth[t - 1, r]) < 0.0
                        v.npwetcost[t, r] = 0
                    else
                        v.npwetcost[t, r] = p.wmbm[r] * ds * v.wetval[t, r] * (1.0 + p.slrprtp[r] + ypcgrowth) / (p.slrprtp[r] + ypcgrowth - p.wvel * ypcgrowth - p.wvpdl * popdensgrowth - p.wvsl * v.wetlandgrowth[t - 1, r])
                    end

                    # NPV of dryland
                    if (1.0 + p.dvydl * incomedensgrowth) < 0.0
                        v.npdrycost[t, r] = 0
                    else
                        v.npdrycost[t, r] = potLandloss * v.dryval[t, r] * (1 + p.slrprtp[r] + ypcgrowth) / (p.slrprtp[r] + ypcgrowth - p.dvydl * incomedensgrowth)
                    end

                    # Calculate protection level
                    v.protlev[t, r] = max(0.0, 1.0 - 0.5 * (v.npprotcost[t, r] + v.npwetcost[t, r]) / v.npdrycost[t, r])

                    if v.protlev[t, r] > 1
                        error("protlevel >1 should not happen")
                    end
                end

                # Calculate actual wetland loss and cost
                v.wetlandloss[t, r] = min(
                    p.wlbm[r] * ds + v.protlev[t, r] * p.wmbm[r] * ds,
                    p.wetmax[r] - v.cumwetlandloss[t - 1, r])

                v.cumwetlandloss[t, r] = v.cumwetlandloss[t - 1, r] + v.wetlandloss[t, r]

                # Calculate wetland growth
                v.wetlandgrowth[t, r] = (p.wetland90[r] - v.cumwetlandloss[t, r]) / (p.wetland90[r] - v.cumwetlandloss[t - 1, r]) - 1.0

                v.wetcost[t, r] = v.wetval[t, r] * v.wetlandloss[t, r]

                v.landloss[t, r] = (1.0 - v.protlev[t, r]) * potLandloss

                v.cumlandloss[t, r] = v.cumlandloss[t - 1, r] + v.landloss[t, r]
                v.drycost[t, r] = v.dryval[t, r] * v.landloss[t, r]

                v.protcost[t, r] = v.protlev[t, r] * p.pc[r] * ds

                if v.landloss[t, r] < 0
                    v.leave[t, r] = 0
                else
                    v.leave[t, r] = p.coastpd[r] * popdens * v.landloss[t, r]
                end

                v.leavecost[t, r] = p.emcst * ypc * v.leave[t, r] / 1000000000
            end

            for destination in d.regions
                enter = 0.0
                for source in d.regions
                    enter += v.leave[t, source] * v.imigrate[source, destination]
                end
                v.enter[t, destination] = enter
            end

            for r in d.regions
                ypc = p.income[t, r] / p.population[t, r] * 1000.0
                v.entercost[t, r] = p.immcst * ypc * v.enter[t, r] / 1000000000
            end
        end
    end
end

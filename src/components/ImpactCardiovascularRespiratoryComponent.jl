@defcomp impactcardiovascularrespiratory begin
    regions = Index()

    basecardvasc = Variable(index=[time,regions])
    baseresp = Variable(index=[time,regions])

    cardheat = Variable(index=[time,regions])
    resp = Variable(index=[time,regions])
    cardcold = Variable(index=[time,regions])

    cardvasc90 = Parameter(index=[regions])
    plus90 = Parameter(index=[regions])
    resp90 = Parameter(index=[regions])
    chplbm = Parameter(index=[regions])
    chmlbm = Parameter(index=[regions])
    chpqbm = Parameter(index=[regions])
    chmqbm = Parameter(index=[regions])
    rlbm = Parameter(index=[regions])
    rqbm = Parameter(index=[regions])
    ccplbm = Parameter(index=[regions])
    ccmlbm = Parameter(index=[regions])
    ccpqbm = Parameter(index=[regions])
    ccmqbm = Parameter(index=[regions])

    plus = Parameter(index=[time,regions])
    temp = Parameter(index=[time,regions])
    urbpop = Parameter(index=[time,regions])
    population = Parameter(index=[time,regions])

    cvlin = Parameter(default=0.025901)
    rlin = Parameter(default=0.001583)
    maxcardvasc = Parameter(default=0.05)

    function run_timestep(p, v, d, t)

        if !is_first(t)
            for r in d.regions
                v.basecardvasc[t, r] = p.cardvasc90[r] + p.cvlin * (p.plus[t, r] - p.plus90[r])
                if v.basecardvasc[t, r] > 1.0
                    v.basecardvasc[t, r] = 1.0
                end

                v.baseresp[t, r] = p.resp90[r] + p.rlin * (p.plus[t, r] - p.plus90[r])
                if v.baseresp[t, r] > 1.0
                    v.baseresp[t, r] = 1.0
                end

                v.cardheat[t, r] = (p.chplbm[r] * p.plus[t, r] + p.chmlbm[r] * (1.0 - p.plus[t, r])) * p.temp[t, r] +
                        (p.chpqbm[r] * p.plus[t, r] + p.chmqbm[r] * (1.0 - p.plus[t, r])) * p.temp[t, r]^2
                v.cardheat[t, r] = v.cardheat[t, r] * p.urbpop[t, r] * p.population[t, r] * 10
                if v.cardheat[t, r] > 1000.0 * p.maxcardvasc * v.basecardvasc[t, r] * p.urbpop[t, r] * p.population[t, r]
                    v.cardheat[t, r] = 1000 * p.maxcardvasc * v.basecardvasc[t, r] * p.urbpop[t, r] * p.population[t, r]
                end
                if v.cardheat[t, r] < 0.0
                    v.cardheat[t, r] = 0
                end

                v.resp[t, r] = p.rlbm[r] * p.temp[t, r] + p.rqbm[r] * p.temp[t, r]^2
                v.resp[t, r] = v.resp[t, r] * p.urbpop[t, r] * p.population[t, r] * 10
                if v.resp[t, r] > 1000 * p.maxcardvasc * v.baseresp[t, r] * p.urbpop[t, r] * p.population[t, r]
                    v.resp[t, r] = 1000 * p.maxcardvasc * v.baseresp[t, r] * p.urbpop[t, r] * p.population[t, r]
                end
                if v.resp[t, r] < 0
                    v.resp[t, r] = 0
                end

                v.cardcold[t, r] = (p.ccplbm[r] * p.plus[t, r] + p.ccmlbm[r] * (1.0 - p.plus[t, r])) * p.temp[t, r] +
                        (p.ccpqbm[r] * p.plus[t, r] + p.ccmqbm[r] * (1.0 - p.plus[t, r])) * p.temp[t, r]^2
                v.cardcold[t, r] = v.cardcold[t, r] * p.population[t, r] * 10
                if v.cardcold[t, r] < -1000 * p.maxcardvasc * v.basecardvasc[t, r] * p.population[t, r]
                    v.cardcold[t, r] = -1000 * p.maxcardvasc * v.basecardvasc[t, r] * p.population[t, r]
                end
                if v.cardcold[t, r] > 0
                    v.cardcold[t, r] = 0
                end
            end
        end
    end
end
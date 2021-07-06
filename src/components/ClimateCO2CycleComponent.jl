@defcomp climateco2cycle begin
    # Anthropogenic CO2 emissions in Mt of C
    mco2 = Parameter(index=[time])

    # Terrestrial biosphere CO2 emissions in Mt of C
    terrestrialco2 = Variable(index=[time])

    # Net CO2 emissions in Mt of C
    globc = Variable(index=[time])

    # Carbon boxes
    cbox = Variable(index=[time,5])

    # Initial carbon box 1
    cbox10 = Parameter()

    # Initial carbon box 2
    cbox20 = Parameter()

    # Initial carbon box 3
    cbox30 = Parameter()

    # Initial carbon box 4
    cbox40 = Parameter()

    # Initial carbon box 5
    cbox50 = Parameter()

    # Carbon decay in box 1
    co2decay1 = Variable()

    # Carbon decay in box 2
    co2decay2 = Variable()

    # Carbon decay in box 3
    co2decay3 = Variable()

    # Carbon decay in box 4
    co2decay4 = Variable()

    # Carbon decay in box 5
    co2decay5 = Variable()

    # Carbon decay in box 1
    lifeco1 = Parameter(default=1.0)

    # Carbon decay in box 2
    lifeco2 = Parameter(default=363.0)

    # Carbon decay in box 3
    lifeco3 = Parameter(default=74.0)

    # Carbon decay in box 4
    lifeco4 = Parameter(default=17.0)

    # Carbon decay in box 5
    lifeco5 = Parameter(default=2.0)

    # Fraction of carbon emission in box 1
    co2frac1 = Parameter(default=0.13)

    # Fraction of carbon emission in box 2
    co2frac2 = Parameter(default=0.2)

    # Fraction of carbon emission in box 3
    co2frac3 = Parameter(default=0.32)

    # Fraction of carbon emission in box 4
    co2frac4 = Parameter(default=0.25)

    # Fraction of carbon emission in box 5
    co2frac5 = Parameter(default=0.1)

    # Atmospheric CO2 concentration
    acco2 = Variable(index=[time])

    # Stock of CO2 in the terrestrial biosphere
    terrco2stock = Variable(index=[time])

    # Temperature
    temp = Parameter(index=[time])

    terrco2sens = Parameter(default=2600.000000000002)
    terrco2stock0 = Parameter(default=1.9e6)

    tempin2010 = Variable()

    function run_timestep(p, v, d, t)

        if is_first(t)
            v.co2decay1 = p.lifeco1
            v.co2decay2 = exp(-1.0 / p.lifeco2)
            v.co2decay3 = exp(-1.0 / p.lifeco3)
            v.co2decay4 = exp(-1.0 / p.lifeco4)
            v.co2decay5 = exp(-1.0 / p.lifeco5)

            v.terrco2stock[t] = p.terrco2stock0

            v.cbox[t,1] = p.cbox10
            v.cbox[t,2] = p.cbox20
            v.cbox[t,3] = p.cbox30
            v.cbox[t,4] = p.cbox40
            v.cbox[t,5] = p.cbox50
            v.acco2[t] = v.cbox[t,1] + v.cbox[t,2] + v.cbox[t,3] + v.cbox[t,4] + v.cbox[t,5]
        else

            if gettime(t) == 2011
                v.tempin2010 = p.temp[TimestepValue(2010)]
            end

            if gettime(t) > 2010
                v.terrestrialco2[t] = (p.temp[t - 1] - v.tempin2010) * p.terrco2sens * v.terrco2stock[t - 1] / p.terrco2stock0
            else
                v.terrestrialco2[t] = 0
            end

            v.terrco2stock[t] = max(v.terrco2stock[t - 1] - v.terrestrialco2[t], 0.0)

            v.globc[t] = p.mco2[t] + v.terrestrialco2[t]

            # Calculate CO2 concentrations
            v.cbox[t,1] = v.cbox[t - 1,1] * v.co2decay1 + 0.000471 * p.co2frac1 * (v.globc[t])
            v.cbox[t,2] = v.cbox[t - 1,2] * v.co2decay2 + 0.000471 * p.co2frac2 * (v.globc[t])
            v.cbox[t,3] = v.cbox[t - 1,3] * v.co2decay3 + 0.000471 * p.co2frac3 * (v.globc[t])
            v.cbox[t,4] = v.cbox[t - 1,4] * v.co2decay4 + 0.000471 * p.co2frac4 * (v.globc[t])
            v.cbox[t,5] = v.cbox[t - 1,5] * v.co2decay5 + 0.000471 * p.co2frac5 * (v.globc[t])

            v.acco2[t] = v.cbox[t,1] + v.cbox[t,2] + v.cbox[t,3] + v.cbox[t,4] + v.cbox[t,5]
        end
    end
end
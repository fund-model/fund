function interact(M, N)
    d = 1.0 + (M * N)^0.75 * 2.01E-5 + (M * N)^1.52 * M * 5.31E-15
    return 0.47 * log(d)
end

@defcomp climateforcing begin
    # Atmospheric CO2 concentration
    acco2 = Parameter(index=[time])

    # Pre-industrial atmospheric CO2 concentration
    co2pre = Parameter()

    # Atmospheric CH4 concentration
    acch4 = Parameter(index=[time])

    # Pre-industrial atmospheric CH4 concentration
    ch4pre = Parameter()

    # Indirect radiative forcing increase for CH4
    ch4ind = Parameter(default=0.4)

    # Atmospheric N2O concentration
    acn2o = Parameter(index=[time])

    # Pre-industrial atmospheric N2O concentration
    n2opre = Parameter()

    # Atmospheric SF6 concentrations
    acsf6 = Parameter(index=[time])

    # Pre-industrial atmospheric SF6 concentration
    sf6pre = Parameter()

    # Radiative forcing from CO2
    rfco2 = Variable(index=[time])

    # Radiative forcing from CH4
    rfch4 = Variable(index=[time])

    # Radiative forcing from N2O
    rfn2o = Variable(index=[time])

    # Radiative forcing from N2O
    rfsf6 = Variable(index=[time])

    # Radiative forcing from SO2
    rfso2 = Parameter(index=[time])

    # Radiative forcing
    radforc = Variable(index=[time])

    # EMF22 radiative forcing
    rfemf22 = Variable(index=[time])

    function run_timestep(p, v, d, t)

        if !is_first(t)
            ch4n2o = interact(p.ch4pre, p.n2opre)

            v.rfco2[t] = 5.35 * log(p.acco2[t] / p.co2pre)

            v.rfch4[t] = 0.036 * (1.0 + p.ch4ind) * (sqrt(p.acch4[t]) - sqrt(p.ch4pre)) - interact(p.acch4[t], p.n2opre) + ch4n2o

            v.rfn2o[t] = 0.12 * (sqrt(p.acn2o[t]) - sqrt(p.n2opre)) - interact(p.ch4pre, p.acn2o[t]) + ch4n2o

            v.rfsf6[t] = 0.00052 * (p.acsf6[t] - p.sf6pre)

            v.radforc[t] = v.rfco2[t] + v.rfch4[t] + v.rfn2o[t] + v.rfsf6[t] + p.rfso2[t]

            v.rfemf22[t] = v.rfco2[t] + v.rfch4[t] + v.rfn2o[t]
        end
    end
end

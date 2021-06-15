"""
Returns one default FUND model and one model with additional emissions of the specified gas in the specified year.
"""
function getmarginalmodels(; gas = :C, emissionyear = 2010, parameters = nothing, yearstorun = 1050)

    # Get default FUND model
    m1 = getfund(nsteps = yearstorun, params = parameters)

    # Get model to add marginal emissions to
    m2 = getfund(nsteps = yearstorun, params = parameters)
    add_comp!(m2, Mimi.adder, :marginalemission, before = :climateco2cycle, first = 1951)
    addem = zeros(yearstorun)
    addem[getindexfromyear(emissionyear)-1:getindexfromyear(emissionyear) + 8] .= 1.0
    set_param!(m2, :marginalemission, :add, addem)

    # Reconnect the appropriate emissions in the marginal model
    if gas == :C
        connect_param!(m2, :marginalemission, :input, :emissions, :mco2)
        connect_param!(m2, :climateco2cycle, :mco2, :marginalemission, :output, repeat([missing], yearstorun + 1))
    elseif gas == :CH4
        connect_param!(m2, :marginalemission, :input, :emissions, :globch4)
        connect_param!(m2, :climatech4cycle, :globch4, :marginalemission, :output, repeat([missing], yearstorun + 1))
    elseif gas == :N2O
        connect_param!(m2, :marginalemission, :input, :emissions, :globn2o)
        connect_param!(m2, :climaten2ocycle, :globn2o, :marginalemission, :output, repeat([missing], yearstorun + 1))
    else
        error("Unknown gas.")
    end

    # Run each model
    run(m1)
    run(m2)

    return m1, m2
end

"""
Returns the social cost per one ton of additional emissions of the specified gas in the specified year. 
Uses the specified eta and prtp for discounting, with the option to use equity weights.
"""
function marginaldamage3(; emissionyear = 2010, parameters = nothing, yearstoaggregate = 1000, gas = :C, useequityweights = false, eta = 1.0, prtp = 0.001)
    yearstorun = min(1050, getindexfromyear(emissionyear) + yearstoaggregate)

    m1, m2 = getmarginalmodels(emissionyear = emissionyear, parameters = parameters, yearstorun = yearstorun, gas = gas)

    damage1 = m1[:impactaggregation, :loss]
    # Take out growth effect effect of run 2 by transforming
    # the damage from run 2 into % of GDP of run 2, and then
    # multiplying that with GDP of run 1
    damage2 = m2[:impactaggregation, :loss] ./ m2[:socioeconomic, :income] .* m1[:socioeconomic, :income]

    # Calculate the marginal damage between run 1 and 2 for each
    # year/region
    marginaldamage = (damage2 .- damage1) / 10000000.0

    ypc = m1[:socioeconomic, :ypc]

    df = zeros(yearstorun + 1, 16)
    if !useequityweights
        for r = 1:16
            x = 1.
            for t = getindexfromyear(emissionyear):yearstorun
                df[t, r] = x
                gr = (ypc[t, r] - ypc[t - 1, r]) / ypc[t - 1,r]
                x = x / (1. + prtp + eta * gr)
            end
        end
    else
        globalypc = m1[:socioeconomic, :globalypc]
        df = Float64[t >= getindexfromyear(emissionyear) ? (globalypc[getindexfromyear(emissionyear)] / ypc[t, r]) ^ eta / (1.0 + prtp) ^ (t - getindexfromyear(emissionyear)) : 0.0 for t = 1:yearstorun + 1, r = 1:16]
    end

    scc = sum(marginaldamage[2:end, :] .* df[2:end, :])

    return scc
end

"""
Returns a matrix of marginal damages per one ton of additional emissions of the specified gas in the specified year.
"""
function getmarginaldamages(; emissionyear=2010, parameters = nothing, yearstoaggregate = 1000, gas = :C) 

    # Calculate number of years to run the models
    yearstorun = min(1050, getindexfromyear(emissionyear) + yearstoaggregate)

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

    # Run both models
    run(m1; ntimesteps = yearstorun)
    run(m2; ntimesteps = yearstorun)

    # Get damages
    damage1 = m1[:impactaggregation, :loss]
    damage2 = m2[:impactaggregation, :loss]

    # Calculate the marginal damage between run 1 and 2 for each year/region
    marginaldamages = (damage2 .- damage1) / 10000000.0

    return marginaldamages
end

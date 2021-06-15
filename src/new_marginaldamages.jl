import Mimi.compinstance

"""
Creates a MarginalModel of FUND with additional emissions in the specified year for the specified gas. 
"""
function create_marginal_FUND_model(; gas = :CO2, emissionyear = 2010, parameters = nothing, yearstorun = 1050)

    # Get default FUND model
    FUND = getfund(nsteps = yearstorun, params = parameters)

    # Build marginal model
    mm = create_marginal_model(FUND)
    m1, m2 = mm.base, mm.modified

    add_marginal_emissions!(m2, emissionyear; gas = gas)

    Mimi.build(m1)
    Mimi.build(m2)
    return mm 
end 

# A component for an emissions pulse to be used in social cost calculations. Computes the `output` vector by adding
#   `add` to `input`. This is similar to the Mimi.adder component, except that it allows missing values to be passed through.
@defcomp emissionspulse begin
    add    = Parameter(index=[time])
    input  = Parameter(index=[time])
    output = Variable(index=[time])

    function run_timestep(p, v, d, t)
        v.output[t] = Mimi.@allow_missing(p.input[t]) + p.add[t]
    end
end

function _gas_normalization(gas::Symbol)
    if gas == :CO2 
        return 12/44    # convert from tons CO2 to tons C
    elseif gas == :CH4
        return 1
    elseif gas == :N2O
        return 28/44    # convert from tons N2O to tons N
    else
        error("Unknown gas :$gas.")
    end
end

"""
Adds an emissionspulse component to `m`, and sets the additional emissions if a year is specified.
The size of the marginal emission pulse can be modified with the `pulse_size` keyword argument, in metric tonnes.
"""
function add_marginal_emissions!(m, year::Union{Int, Nothing} = nothing; gas::Symbol = :CO2, pulse_size::Float64 = 1e7)

    # Add additional emissions to m
    add_comp!(m, emissionspulse, before = :climateco2cycle)
    nyears = length(Mimi.time_labels(m))
    addem = zeros(nyears) 
    if year !== nothing 
        # pulse is spread over ten years, and emissions components is in Mt so divide by 1e7, and convert from CO2 to C if gas==:CO2 because emissions component is in MtC
        addem[getindexfromyear(year):getindexfromyear(year) + 9] .= pulse_size / 1e7 * _gas_normalization(gas)
    end
    set_param!(m, :emissionspulse, :add, addem)

    # Reconnect the appropriate emissions in m
    if gas == :CO2
        connect_param!(m, :emissionspulse, :input, :emissions, :mco2)
        connect_param!(m, :climateco2cycle, :mco2, :emissionspulse, :output)
    elseif gas == :CH4
        connect_param!(m, :emissionspulse, :input, :emissions, :globch4)
        connect_param!(m, :climatech4cycle, :globch4, :emissionspulse, :output)
    elseif gas == :N2O
        connect_param!(m, :emissionspulse, :input, :emissions, :globn2o)
        connect_param!(m, :climaten2ocycle, :globn2o, :emissionspulse, :output)
    else
        error("Unknown gas: $gas")
    end

end 

"""
Helper function to set the marginal emissions in the specified year.
"""
function perturb_marginal_emissions!(m::Model, year; comp_name::Symbol = :emissionspulse, pulse_size::Float64 = 1e7, gas::Symbol = :CO2)

    ci = compinstance(m, comp_name)
    emissions = Mimi.get_param_value(ci, :add)

    nyears = length(Mimi.dimension(m, :time))
    new_em = zeros(nyears)
    new_em[getindexfromyear(year):getindexfromyear(year) + 9] .= pulse_size / 1e7 * _gas_normalization(gas)
    emissions[:] = new_em

end

"""
Returns the social cost per one ton of additional emissions of the specified gas in the specified year. 
Uses the specified eta and prtp for discounting, with the option to use equity weights.
"""
function get_social_cost(; emissionyear = 2010, parameters = nothing, yearstoaggregate = 1000, gas = :CO2, useequityweights = false, eta = 1.0, prtp = 0.001)

    # Get marginal model
    yearstorun = min(1050, getindexfromyear(emissionyear) + yearstoaggregate)
    mm = create_marginal_FUND_model(emissionyear = emissionyear, parameters = parameters, yearstorun = yearstorun, gas = gas)
    run(mm)
    m1, m2 = mm.base, mm.modified

    damage1 = m1[:impactaggregation, :loss]
    # Take out growth effect effect of run 2 by transforming the damage from run 2 into % of GDP of run 2, and then multiplying that with GDP of run 1
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

"""
Returns a matrix of marginal damages per one ton of additional emissions of the specified gas in the specified year.
"""
function getmarginaldamages(; emissionyear=2010, parameters = nothing, yearstoaggregate = 1000, gas = :CO2) 

    # Get marginal model
    yearstorun = min(1050, getindexfromyear(emissionyear) + yearstoaggregate)
    mm = create_marginal_FUND_model(emissionyear = emissionyear, parameters = parameters, yearstorun = yearstorun, gas = gas)
    run(mm)

    # Get damages
    marginaldamages = mm[:impactaggregation, :loss] / 10000000.0
    return marginaldamages
end
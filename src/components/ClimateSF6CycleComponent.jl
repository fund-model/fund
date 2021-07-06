@defcomp climatesf6cycle begin
    # Global SF6 emissions in kt of SF6
    globsf6 = Parameter(index=[time])

    # Atmospheric SF6 concentrations
    acsf6 = Variable(index=[time])

    # SF6 pre industrial
    sf6pre = Parameter()

    # SF6 decay
    sf6decay = Variable()

    lifesf6 = Parameter(default=3200.0)

    function run_timestep(p, v, d, t)
        
        if is_first(t)
            v.sf6decay = 1.0 / p.lifesf6

            v.acsf6[t] = p.sf6pre
        else
            # Calculate SF6 concentrations
            v.acsf6[t] = p.sf6pre + (v.acsf6[t - 1] - p.sf6pre) * (1 - v.sf6decay) + p.globsf6[t] / 25.1

            if v.acsf6[t] < 0
                error("sf6 atmospheric concentration out of range")
            end
        end
    end
end

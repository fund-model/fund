@defcomp climatech4cycle begin
    # Global CH4 emissions in Mt of CH4
    globch4 = Parameter(index=[time])

    # Atmospheric CH4 concentration
    acch4 = Variable(index=[time])

    # CH4 decay
    ch4decay = Variable()

    lifech4 = Parameter(default=12.0)

    #  CH4 pre industrial
    ch4pre = Parameter()

    # Initial acch4 value
    acch4_0 = Parameter(default=1222.0)

    function run_timestep(p, v, d, t)
        
        if is_first(t)
            v.ch4decay = 1.0 / p.lifech4

            v.acch4[t] = p.acch4_0
        else
            # Calculate CH4 concentrations
            v.acch4[t] = v.acch4[t - 1] + 0.3597 * p.globch4[t] - v.ch4decay * (v.acch4[t - 1] - p.ch4pre)

            if v.acch4[t] < 0
                error("ch4 atmospheric concentration out of range in $t")
            end
        end
    end
end

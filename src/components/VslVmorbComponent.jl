@defcomp vslvmorb begin
    regions = Index()

    vsl = Variable(index=[time,regions])
    vmorb = Variable(index=[time,regions])

    population = Parameter(index=[time,regions])
    income = Parameter(index=[time,regions])

    vslbm       = Parameter(default=4.99252262888626e6)
    vslel       = Parameter(default=1)
    vmorbbm     = Parameter(default=19970.090515545)
    vmorbel     = Parameter(default=1)
    vslypc0     = Parameter(default=24962.6131444313)
    vmorbypc0   = Parameter(default=24962.6131444313)

    function run_timestep(p, v, d, t)

        if !is_first(t)
            for r in d.regions
                ypc = p.income[t, r] / p.population[t, r] * 1000.0

                v.vsl[t, r] = p.vslbm * (ypc / p.vslypc0)^p.vslel

                v.vmorb[t, r] = p.vmorbbm * (ypc / p.vmorbypc0)^p.vmorbel
            end
        end
    end
end
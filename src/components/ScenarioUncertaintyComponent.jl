@defcomp scenariouncertainty begin
    regions = Index()

    pgrowth = Variable(index=[time,regions])
    ypcgrowth = Variable(index=[time,regions])
    aeei = Variable(index=[time,regions])
    acei = Variable(index=[time,regions])
    forestemm = Variable(index=[time,regions])

    timeofuncertaintystart = Parameter{Int}(default=2000)

    scenpgrowth     = Parameter(index=[time,regions])
    scenypcgrowth   = Parameter(index=[time,regions])
    scenaeei        = Parameter(index=[time,regions])
    scenacei        = Parameter(index=[time,regions])
    scenforestemm   = Parameter(index=[time,regions])

    ecgradd     = Parameter(index=[regions])
    pgadd       = Parameter(index=[regions])
    aeeiadd     = Parameter(index=[regions])
    aceiadd     = Parameter(index=[regions])
    foremadd    = Parameter(index=[regions])

    function run_timestep(p, v, d, t)
        yearsFromUncertaintyStart = gettime(t) - p.timeofuncertaintystart
        sdTimeFactor = (yearsFromUncertaintyStart / 50.0) / (1.0 + (yearsFromUncertaintyStart / 50.0))

        for r in d.regions
            v.ypcgrowth[t, r] = p.scenypcgrowth[t, r] + (gettime(t) >= p.timeofuncertaintystart ? p.ecgradd[r] * sdTimeFactor : 0.0)
            v.pgrowth[t, r] = p.scenpgrowth[t, r] + (gettime(t) >= p.timeofuncertaintystart ? p.pgadd[r] * sdTimeFactor : 0.0)
            v.aeei[t, r] = p.scenaeei[t, r] + (gettime(t) >= p.timeofuncertaintystart ? p.aeeiadd[r] * sdTimeFactor : 0.0)
            v.acei[t, r] = p.scenacei[t, r] + (gettime(t) >= p.timeofuncertaintystart ? p.aceiadd[r] * sdTimeFactor : 0.0)
            v.forestemm[t, r] = p.scenforestemm[t, r] + (gettime(t) >= p.timeofuncertaintystart ? p.foremadd[r] * sdTimeFactor : 0.0)
        end
    end
end
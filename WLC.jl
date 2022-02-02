using Roots

"""
Return y'(x) using forward numerical differentiation.
"""
function num_differentiate(xs,ys)
    ders=zeros(length(xs)-1)
    for i in 1:length(xs)-1
        ders[i]=(ys[i+1]-ys[i])/(xs[i+1]-xs[i])
    end
    return ders
end

"""
Function where root x is relative extension of
worm-like chain in l>P regime with
    F = force
    P = persistence length

Bouchiat C, Wang MD, Allemand J, Strick T, Block SM, Croquette V.
Estimating the persistence length of a worm-like chain molecule from force-extension measurements.
Biophys J. 1999;76:409-413.
doi:10.1016/s0006-3495(99)77207-3
"""
function g_long(x, F, P, kT)
    a_two=-0.5164228
    a_three=-2.737418
    a_four=16.07497
    a_five=-38.87607
    a_six=39.49944
    a_seven=-14.17718

    return F*P/kT -(1-x)^(-2)/4. + 1. /4. -x - a_two*x^2-a_three*x^3-a_four*x^4-a_five*x^5-a_six*x^6-a_seven*x^7
end

"""
Function where root x is relative extension of
worm-like chain in other regime with l<P
    F = force
    P = persistence length
    kT = k_B T

Broedersz C.P., MacKintosh F.C.
Modeling semiflexible polymer networks
Rev. Mod. Phys. 2014;86:995-1036
doi:10.1103/RevModPhys.86.995
"""
function g_short(x, F, P, l, kT)
    ϕ=F*l^2/(pi^2*P*kT)
    return ϕ-9/pi^2*((1-x)^(-2)-1-x/3)
end

"""
Finds the relative extension of a chain with
    F = force
    P = persistence length
    kT = k_B T
    Marko_Siggia : if true, use Marko-Siggia, else
"""
function wlc_ext(L,P,F,kT=1, Marko_Siggia=true)
    if Marko_Siggia
        return L*find_zero(x->g_long(x,F,P,kT),(0,1),Bisection())
    else
        return L*find_zero(x->g_short(x,F,P,L,kT),(0,1),Bisection())
    end

end

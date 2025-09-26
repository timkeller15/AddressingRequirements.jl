function reflection_coefficient_simple(N::Int64, atom::Atom, cav::Cavity)

    r  = 1. - 2*cav.κr/(cav.κ + N*atom.coupling^2/atom.γ) 

    return r
end

function fidelity_local_ana(coupling::Float64, cav::Cavity)

    atom = Atom(coupling = coupling)
    r0 = reflection_coefficient_simple(0,atom,cav)
    r1 = reflection_coefficient_simple(1,atom,cav) 
    r2 = reflection_coefficient_simple(2,atom,cav) 

    F = 0.25*abs2(r0 - 2*r1 - r2)/(abs2(r0) + 2*abs2(r1) + abs2(r2))

    return F
end

function fidelity_local_ana_full(cavity::Cavity, sys::System)

    N = sys.N

    selection = [i+1 for i = 0:(2^N-1) if mod(floor(i/2^(sys.targets[1]-1)),2)==0 && mod(floor(i/2^(sys.targets[2]-1)),2)==0]
    Fij = ones(2^N)
    Fij[selection] .= -1
    Fij = Fij*Fij'

    G = G_matrix(cavity,sys)
    F = real(sum(G.*Fij)/(2^N*tr(G)))

    return F
end

function fidelity_remote_ana(coupling::Float64, cav::Cavity)

    atom = Atom(coupling = coupling)
    r0 = reflection_coefficient_simple(0,atom,cav) 
    r1 = reflection_coefficient_simple(1,atom,cav) 

    F_H = 0.25*abs2(r0^2 + 4*r0 + 2*r1*r0 - r1^2 - 2)/(abs2(r0^2 + 2*r0 - 1) + 2*abs2(r1*r0 + r1 + r0 - 1) + abs2(r1^2 + 2*r1 - 1))
    F_V = 0.25*abs2(r0^2 - 2*r0 + r1^2 + 2*r1 + 2)/(abs2(r0^2 + 1) + abs2(r1*r0 + r1 - r0 + 1) + abs2(r1*r0 + r0 - r1 + 1) + abs2(r1^2 + 1))

    F = 0.5*(F_H + F_V)

    return F
end
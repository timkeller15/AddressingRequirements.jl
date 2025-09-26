function collective_coupling_strength(cavity::Cavity, sys::System) 
    g_arr = []
    for atom in cavity.atoms
        push!(g_arr,atom.coupling^2/(1im*(sys.ωp - atom.ωa) - atom.γ))
    end

    g = sum(g_arr)

    return g
end

function reflection_coefficient(cavity::Cavity, sys::System) 

    g = collective_coupling_strength(cavity,sys)
    r = 1 + 2*cavity.κr/(1im*(sys.ωp - cavity.ωc) - cavity.κ + g)

    return r
end

function G_matrix(cavity::Cavity, sys::System) 
    
    N = sys.N
    p = Complex.(zeros(2^N))

    for i = 1:2^N

        cav = deepcopy(cavity)
        configuration = digits(i-1,base=2,pad=N)

        for atom in cav.atoms
            atom.ωa = atom.ωa + atom.ωq*(1 - configuration[atom.ID])
            
            # adjust coupling strength of non-coupling state to physical Cesium setup by factor of sqrt(3/7) stemming from difference f = -sqrt(5/18) for |1> and f = -sqrt(5/42) for |0> state
            if sys.adjust_coupling
                coupling_adjusted = atom.coupling*(configuration[atom.ID] + (1 - configuration[atom.ID])*sqrt(3/7))
                set_coupling(atom,coupling_adjusted)
            end
        end

        p[i] = reflection_coefficient(cav,sys)
    end

    if sys.mode == :remote
        p = kron(p,[1,0]) + kron(ones(2^N),[0,1])
    end

    G = p*p'

    return G
end

function success_probability(sys::System)
    if sys.mode == :local
        return success_probability_local(sys)
    elseif sys.mode == :remote
        return success_probability_remote(sys)
    else
        return nothing
    end
end

function success_probability_local(sys::System; cavity::Cavity = sys.cavities[1])
    G = G_matrix(cavity,sys) 
    return real(tr(G)/2^sys.N)
end

function success_probability_remote(sys::System)
    N = sys.N

    Hadamard = [1 1; 1 -1]/sqrt(2)
    Hadamard = kron(collect(I(2^N)),Hadamard)

    G1 = G_matrix(sys.cavities[1],sys)   
    G2 = G_matrix(sys.cavities[2],sys) 

    GG = Hadamard*((Hadamard*G1*Hadamard').*G2)*Hadamard'
    
    return real(tr(GG)/2^(N+1))
end

density_matrix(state::Vector{Float64}) = state*state'

function choi(input::Matrix{T}) where {T<:Number}

    N = Int(log2(length(input))/2)

    return kron(collect(I(2^N)),input)
end

function choi_state(N::Int64)
        
    Φ⁺ = dropdims(reshape(collect(I(2^N)),(2^(2*N),1)),dims=2)/sqrt(2^N)
    # alternatively: Φ⁺ = sum([(input = zeros(2^N); input[i] = 1; kron(input,input)) for i=1:2^N])/sqrt(2^N)

    return Φ⁺
end

function fidelity(input::Matrix{ComplexF64}, target::Vector{Float64}; choi::Bool = false)

    fidelity = (target'*input*target)/tr(input)
    fidelity = real(fidelity)

    if choi
        N = Int(log2(length(target))/2)
        fidelity = (2^N*fidelity + 1)/(2^N + 1)
    end

    return fidelity
end
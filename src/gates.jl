function CZ_unitary(input::Union{Vector{Float64},Nothing}, N::Int64, control::Int64, target::Int64; mode::Symbol)

    if typeof(input) == Nothing
        input = 1
    end

    configurations = digits.(0:2^N-1,base=2,pad=N)

    if mode == :default
        selection = [i for i = 1:length(configurations) if configurations[i][control]==1 && configurations[i][target]==1]
    elseif mode == :phase_flip
        selection = [i for i = 1:length(configurations) if configurations[i][control]==0 && configurations[i][target]==0]
    else
        selection = []
    end

    U_CZ = ones(2^N)
    U_CZ[selection] .= -1
    U_CZ = diagm(U_CZ)

    return U_CZ*input
end

function CZ_unitary(input::Union{Vector{Float64},Nothing}, sys::System)

    if sys.mode == :local
        mode = :phase_flip
    elseif sys.mode == :remote
        mode = :default
    else
        mode = :none
    end

    return CZ_unitary(input, sys.N, sys.targets[1], sys.targets[2]; mode = mode)
end

function local_gate(input::Union{Matrix{T},Nothing}, sys::System; cavity::Cavity = sys.cavities[1], choi::Bool = false, error_rates::Bool = false) where {T<:Number}

    G = G_matrix(cavity,sys)    

    if choi
        N = sys.N
        Φ⁺ = choi_state(N)
        rho = kron(collect(ones(2^N,2^N)),G).*(Φ⁺*Φ⁺')
    else
        rho = G.*input
    end

    if error_rates && !choi 
        U = CZ_unitary(nothing,sys)
        rho = U*rho*U'
    end

    return rho
end

function remote_gate(input::Union{Matrix{T},Nothing}, sys::System; choi::Bool = false, error_rates::Bool = false) where {T<:Number}
   
    N = sys.N

    photon = [1,1]/sqrt(2)
    
    Hadamard = [1 1; 1 -1]/sqrt(2)
    Hadamard = kron(collect(I(2^N)),Hadamard)

    σz = [1 0; 0 -1]
    # Z-Gate on Qubit 1 (little endian!)
    σz = kron(collect(I(2^(N-sys.targets[1]))),σz,collect(I(2^(sys.targets[1]-1)))) 
    
    G1 = G_matrix(sys.cavities[1],sys)   
    G2 = G_matrix(sys.cavities[2],sys)   

    horz = kron(collect(I(2^N)),[1,0]) 
    vert = kron(collect(I(2^N)),[0,1]) 

    if choi
        Φ⁺ = choi_state(N)
        Φ⁺ = kron(Φ⁺,photon)

        Hadamard = kron(collect(I(2^N)),Hadamard)
        σz = kron(collect(I(2^N)),σz) 
        horz = kron(collect(I(2^N)),horz)
        vert = kron(collect(I(2^N)),vert)
        
        rho = kron(collect(ones(2^N,2^N)),G1).*(Φ⁺*Φ⁺')
        rho = Hadamard*rho*Hadamard'
        rho = kron(collect(ones(2^N,2^N)),G2).*rho
        rho = Hadamard*rho*Hadamard'

    else
        rho = kron(input,density_matrix(photon))
        
        rho = G1.*rho
        rho = Hadamard*rho*Hadamard'
        rho = G2.*rho
        rho = Hadamard*rho*Hadamard'
    end

    # Z-Gate on Qubit 1 (little endian!)
    rho_H = (horz')*rho*horz

    rho_V = (vert')*rho*vert
    rho_V = σz'*rho_V*σz

    # tracing out the photon
    rho = rho_H + rho_V

    if error_rates && !choi 
        U = CZ_unitary(nothing,sys)
        rho = U*rho*U'
    end

    return rho
end
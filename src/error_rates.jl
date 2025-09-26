function pauli_basis(N::Int64)
    σ0 = [1 0; 0 1]; 
    σx = [0 1; 1 0];
    σy = [0 -1im; 1im 0]; 
    σz = [1 0; 0 -1]
    
    basis = [σ0,σx,σy,σz]

    if N > 1
        return [kron(σi,σj) for σi in basis for σj in pauli_basis(N-1)]
    else
        return basis
    end
end

function QGT_matrix(N::Int64)
    QGT = [1 1 1 1; 1 1 -1 -1; 1 -1 1 -1; 1 -1 -1 1]/4

    if N > 1
        return kron(QGT,QGT_matrix(N-1))
    else
        return QGT 
    end
end

function pauli_twirling(channel::Function, input::Matrix{T}, sys::System; mode::Symbol = :average, kwargs...) where {T<:Number}  

    basis = pauli_basis(sys.N)

    if mode == :average
        χ = Complex.(zeros(size(input)))
        for σn in basis
            χ += σn*channel(σn*input*σn,sys; kwargs...)*σn
        end
        χ /= length(basis)
    end

    if mode == :random
        σn = rand(basis)
        χ = σn*channel(σn*input*σn,sys)*σn
    end
    
    return χ
end

function pauli_transfer_matrix(channel::Function, sys::System; kwargs...)
    # calculate the diagonal of the Pauli transfer matrix
    N = sys.N
    M = zeros(4^N)
    
    for (i,σn) in enumerate(pauli_basis(N))
        input = σn 
        χ = pauli_twirling(channel,input,sys; kwargs...)
        M[i] = real(tr(σn*χ)/2^N)
    end
    
    return M
end

function gate_error_rates(channel::Function, sys::System; trace_out::Bool = false)
    error_rates = QGT_matrix(sys.N)*pauli_transfer_matrix(channel,sys; error_rates = true)

    if trace_out
       error_rates = sum(reshape(error_rates,(4^2,4^(sys.N-2))),dims=2)
    end

    return error_rates
end

function bias(probabilities)
    
    η = NaN
    
    if length(probabilities) == 4
        η = p[4]/(p[2] + p[3]) # for single qubit: pz/(px + py)
    elseif length(probabilities) == 16
        dephasing = [4,13,16] # IZ, ZI, ZZ errors
        non_dephasing = [2,3,5,6,7,8,9,10,11,12,14,15] # rest 
        p_dephasing = sum(probs[dephasing])
        p_non_dephasing = sum(probs[non_dephasing])
        
        η = p_dephasing/p_non_dephasing 
    end
    
    return η
end
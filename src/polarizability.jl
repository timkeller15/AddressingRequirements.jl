struct Level
    name::String
    n::Int64
    l::Int64
    s::Float64
    J::Float64
    mJ::Union{Float64,Missing}
    F::Union{Float64,Missing}
    mF::Union{Float64,Missing}
    type::Symbol

    function Level(name::String, n::Int64, l::Union{Int64,String}, J::Float64; s::Float64 = 0., mJ::Union{Float64,Missing} = missing, F::Union{Float64,Missing} = missing, mF::Union{Float64,Missing} = missing)
        mapping = Dict("S"=>0, "P"=>1, "D"=>2, "F"=>3, "G"=>4, "H"=>5)
        if typeof(l) == String
            l = get(mapping, uppercase(l), 0)
        end
        s = J - l
    
        if typeof(mJ) == Float64
            type = :fine
        elseif typeof(F) == Float64 && typeof(mF) == Float64
            type = :hyperfine
        else
            type = :atomic
        end

        return new(name,n,l,s,J,mJ,F,mF,type)
    end
end

function Level(input::String; mJ::Union{Real,Missing} = missing, F::Union{Real,Missing} = missing, mF::Union{Real,Missing} = missing)
    
    types = [isdigit(char) for char in input]
    ind1 = findfirst(!identity,types)
    ind2 = findlast(!identity,types)

    n = parse(Int64,input[1:ind1-1])
    l = string(input[ind1])
    J = parse(Float64,input[ind1+1:ind2-1])/parse(Float64,input[ind2+1:end])
    
    return Level(input,n,l,J; mJ = float(mJ), F = float(F), mF = float(mF))
end

function Base.show(io::IO, level::Level) 
    println(io, titlecase(String(getfield(level,:type))), " Level ", getfield(level,:name)) 
end

struct Transition
    initial::Level
    final::Level
    λ::Float64
    ω::Float64
    dipole::Float64
    γ::Float64

    function Transition(initial::Level, final::Level, λ::Float64, dipole::Float64, γ::Union{Float64,Missing}) 
        
        ω = 2*π*c/λ

        if typeof(γ) == Missing
            cj = (2*final.J + 1)/(2*initial.J + 1)
            γ = 0.5*cj*ω^3*abs2(e*a0*dipole)/(3*π*ϵ0*ħ*c^3)  # https://steck.us/alkalidata/
            γ /= 2*π*1e6
        end

        return new(initial,final,λ,ω,dipole,γ)
    end
end

function Base.show(io::IO, transition::Transition) 
    println(io, "Transition ", getfield(getfield(transition,:initial),:name), " --> ", getfield(getfield(transition,:final),:name), " | λ = ", @sprintf("%.1f nm",1e9*getfield(transition,:λ)), " | Dipole Moment ", @sprintf("%2.1f ea0",getfield(transition,:dipole)), " | γ = ", @sprintf("%2.1f MHz",getfield(transition,:dipole))) 
end

struct LevelStructure
    atom::Symbol
    I::Union{Float64,Missing}
    nr_levels::Int64
    nr_transitions::Int64
    levels::Vector{Level}
    transitions::Vector{Transition}
    
    function LevelStructure(atom::Symbol)
        
        levels = Level[]
        transitions = Transition[]
        
        nuclear_spin = Dict(:cesium => 7/2)
        I = get(nuclear_spin,atom,missing)

        data = readdlm("data/" * lowercase(String(atom)) * "_levelstructure.dat", comments = true)

        for i = 1:size(data)[1]
            initial = Level(String(data[i,1]))
            final = Level(String(data[i,2]))

            λ = float(data[i,3])*1e-9
            dipole = float(data[i,4])

            if size(data)[2] > 4 && typeof(data[i,5]) == Float64
                γ = data[i,5]
            else
                γ = missing
            end

            transition = Transition(initial,final,λ,dipole,γ)
            push!(transitions,transition)
            push!(levels,initial,final)
        end
        
        levels = sort(unique(levels), by = lvl -> (lvl.n,lvl.J))
            
        return new(atom,I,length(levels),length(transitions),levels,transitions)
    end
end

function Base.show(io::IO, LS::LevelStructure) 
    println(io, uppercasefirst(String(getfield(LS,:atom))), " Level Structure | I = ", getfield(LS,:I), " | ", getfield(LS,:nr_levels), " Levels | ", getfield(LS,:nr_transitions), " Transitions")
end

function convert_alpha(α::Float64, to::Symbol)
    
    factor = μe*(e/(4*π*ϵ0*ħ*a0))^2 # from Eur. Phys. J. D 67, 92 (2013).
    
    if to == :AU
        exp = 1
    elseif to == :SI
        exp = -1
    else 
        exp = 0
    end
    
    return α*factor^exp 
end

function reduced_alpha(K::Int64, level::Level, LS::LevelStructure, ω::Float64)
    
    arr = zeros(length(LS.transitions)) 
        
    for (i,transition) in enumerate(LS.transitions)
        if transition.final.name == level.name
            J = transition.final.J
            JP = transition.initial.J
            dipole = e*a0*transition.dipole
            γ = 2*π*1e6*transition.γ
            ω1 = (transition.ω - ω - 0.5*1im*γ)
            ω2 = (transition.ω + ω + 0.5*1im*γ)
            sign = Complex(-1)^(K + J + 1 + JP)
            arr[i] = sign*wigner6j(Float64, 1, K, 1, J, JP, J)*abs2(dipole)*real(1/ω1 + (-1)^K/ω2)
        end
    end
    
    α = sqrt(2*K + 1)*sum(arr) 
    
    return α/ħ
end

function polarizability(level::Level, LS::LevelStructure, ω::Float64) 
    
    θp = 0. # angle between quantization axis and polarization
    θk = 0. # angle between quantization axis and propagation direction
    χ = 0. # polarization ellipticity parameter 

    I = LS.I
    J = level.J    

    αS = reduced_alpha(0,level,LS,ω)/sqrt(3*(2*J + 1))
        
    if level.type == :fine
        mJ = level.mJ

        αI = -sqrt(2*J/((J + 1)*(2*J + 1)))*reduced_alpha(1,level,LS,ω) # can be left out when only considering linear polarization
        αT = -sqrt(2*J*(2*J - 1)/(3*(J + 1)*(2*J + 1)*(2*J + 3)))*reduced_alpha(2,level,LS,ω)

        α = αS + sin(2*χ)*cos(θk)*0.5*(mJ/J)*αI + ((3*cos(θp)^2 - 1)/2)*(3*mJ^2 - J*(J + 1))/(J*(2*J - 1))*αT 
    elseif level.type == :hyperfine
        F = level.F
        mF = level.mF

        αI = (-1)^(J + I + F)*sqrt(2*F*(2*F + 1)/(F + 1))*wigner6j(Float64, F, 1, F, J, I, J)*reduced_alpha(1,level,LS,ω) # can be left out when only considering linear polarization
        αT = (-1)^(J + I + F + 1)*sqrt(2*F*(2*F - 1)*(2*F + 1)/(3*(F + 1)*(2*F + 3)))*wigner6j(Float64, F, 2, F, J, I, J)*reduced_alpha(2,level,LS,ω)

        α = αS + sin(2*χ)*cos(θk)*0.5*(mF/F)*αI + ((3*cos(θp)^2 - 1)/2)*(3*mF^2 - F*(F + 1))/(F*(2*F - 1))*αT 
    else
        α = αS
    end

    α_core = convert_alpha(15.84 - 0.67,:SI) # from Phys. Rev. A 104, 052813 (2021)
    α += α_core

    return α
end
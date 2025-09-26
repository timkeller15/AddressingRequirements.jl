mutable struct Atom
    ID::Union{Int64, Missing}
    ωq::Float64
    γ::Float64
    coupling::Float64
    ωa::Float64 
    target::Bool

    Atom(; ID::Union{Int64, Missing} = missing, coupling::Float64 = 0., ωa::Float64 = 0., ωq::Float64 = 9.1926e3, γ::Float64 = 2.6, target::Bool = false) = new(ID,ωq,γ,coupling,ωa,target)
end

function Base.show(io::IO, atom::Atom) 
    if typeof(getfield(atom,:ID)) == Missing
        name = "Atom"
    else
        name = "Atom " * @sprintf("%d",getfield(atom,:ID))
    end
    println(io, name, " | Target: ", getfield(atom,:target), " | ωq = ", getfield(atom,:ωq), " | γ = ", getfield(atom,:γ), " | g = ", getfield(atom,:coupling), " | ωa = ", getfield(atom,:ωa))
end

mutable struct Cavity
    ID::Union{Int64, Missing}
    ωc::Float64
    κr::Float64
    κm::Float64
    κt::Float64
    κ::Float64
    atoms::Vector{Atom}

    Cavity(; ID::Union{Int64, Missing} = missing, ωc::Float64 = 0., κr::Float64 = 2.5, κm::Float64 = 0.1, κt::Float64 = 0.1, atoms::Vector{Atom} = Vector{Atom}()) = new(ID,ωc,κr,κm,κt,sum(κr + κm + κt),atoms)
end

function Base.show(io::IO, cav::Cavity) 

    if typeof(getfield(cav,:ID)) == Missing
        name = "Cavity"
    else
        name = "Cavity " * @sprintf("%d",getfield(cav,:ID))
    end
    
    layout = ""
    for atom in cav.atoms
        if atom.target
            layout *= "x"
        else
            layout *= "o"
        end
    end
    layout = "[" * layout * "]"
    decay_rates = "(κr,κm,κt) = (" * @sprintf("%2.2f",getfield(cav,:κr)) * "," * @sprintf("%2.2f",getfield(cav,:κm)) * "," * @sprintf("%2.2f",getfield(cav,:κm)) * ")"
    println(io, name, " ", layout," | κr/κ = " , @sprintf("%2.4f",getfield(cav,:κr)/getfield(cav,:κ)) , " | ", decay_rates, " | ωc = ", getfield(cav,:ωc))
end

mutable struct System
    N::Int64
    mode::Symbol
    targets::Vector{Int64}
    cavities::Vector{Cavity}
    ωp::Float64
    adjust_coupling::Bool
end

function System(; N::Int64 = 2, mode::Symbol = :local, targets::Vector{Int64} = [1,2], ωp::Float64 = 0., coupling::Float64 = 0., κr::Float64 = 2.5, κm::Float64 = 0.1, κt::Float64 = 0.1, ωq::Float64 = 9.1926e3, adjust_coupling::Bool = false)
   
    if mode == :local
        cavities = [Cavity(ID = 1; κr, κm, κt)]
    elseif mode == :remote
        cavities = [Cavity(ID = 1; κr, κm, κt), Cavity(ID = 2; κr, κm, κt)]
    else
        cavities = []
    end

    sys = System(N,mode,targets,cavities,ωp,adjust_coupling)
    distribute_atoms(sys,ωq)
    set_coupling(sys,coupling)
    set_cavity_decay(sys,κr)

    return sys
end

function Base.show(io::IO, sys::System)

        layout = ""

        for cavity in sys.cavities
            layout *= "["
            for atom in cavity.atoms
                if atom.target
                    layout *= "x"
                else
                    layout *= "o"
                end
            end
            layout *= "]--"
        end
        if getfield(sys,:mode) == :remote
            layout = layout[1:end-2]
        end
  
        println(io, uppercasefirst(String(getfield(sys,:mode))), " Gate | N = ", getfield(sys,:N), " | ", layout, " | ωp = ", getfield(sys,:ωp))
end

function distribute_atoms(sys::System, ωq::Float64)
    
    for i = 1:sys.N
        ind = mod(i-1,length(sys.cavities)) + 1
        atom = Atom(ID = i, ωq = ωq)
        push!(sys.cavities[ind].atoms,atom)

        if i in sys.targets
            atom.target = true
        end
    end
end

function get_targets(sys::System)
    targets = []
    for cavity in sys.cavities
        for atom in cavity.atoms
            if atom.target
                push!(targets,atom)
            end
        end
    end
    return targets
end

function set_coupling(atom::Atom, coupling::Float64)
    atom.coupling = coupling
    return nothing 
end

function set_coupling(cavity::Cavity, coupling::Float64)
    for atom in cavity.atoms
        set_coupling(atom,coupling)
    end
    return nothing 
end

function set_coupling(sys::System, coupling::Float64)
    for cavity in sys.cavities
        set_coupling(cavity,coupling)
    end
    return nothing 
end

function set_detuning(atom::Atom, detuning::Float64)
    if !atom.target
        atom.ωa = detuning
    end
    return nothing 
end

function set_detuning(cavity::Cavity, detuning::Float64)
    for atom in cavity.atoms
        set_detuning(atom,detuning)
    end
    return nothing 
end

function set_detuning(sys::System, detuning::Float64)
    for cavity in sys.cavities
        set_detuning(cavity,detuning)
    end
    return nothing 
end

function set_cavity_decay(cavity::Cavity, cavity_decay::Float64)
    cavity.κr = cavity_decay
    cavity.κ = sum(cavity.κr + cavity.κm + cavity.κt)
    return nothing 
end

function set_cavity_decay(sys::System, cavity_decay::Float64)
    for cavity in sys.cavities
        set_cavity_decay(cavity,cavity_decay)
    end
    return nothing 
end
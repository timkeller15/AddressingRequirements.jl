module AddressingRequirements

    using LinearAlgebra, SpecialFunctions, NonlinearSolve, WignerSymbols
    using Printf, DelimitedFiles

    include("structs.jl")
    export Atom, Cavity, System
    export distribute_atoms, get_targets
    export set_coupling, set_detuning, set_cavity_decay

    include("gates.jl")
    export CZ_unitary, local_gate, remote_gate

    include("error_rates.jl")
    export pauli_basis, QGT_matrix, pauli_twirling, pauli_transfer_matrix, gate_error_rates, bias

    include("analytics.jl")
    export reflection_coefficient_simple, fidelity_local_ana, fidelity_local_ana_full, fidelity_remote_ana

    include("observables.jl")
    export collective_coupling_strength, reflection_coefficient, G_matrix
    export success_probability, success_probability_local, success_probability_remote
    export density_matrix, choi, choi_state
    export fidelity

    include("nanofiber_modes.jl")
    export fiber_parameters, fiber_mode, fiber_mode_linear, fiber_coupling

    include("polarizability.jl")
    export LevelStructure, Level, polarizability, convert_alpha

    include("helper.jl")
    export c, ħ, ϵ0, save_data
end 
# AdressingRequirements.jl

This repository contains a Julia module of the code and data necessary for recreating the findings in the publication

*Addressing requirements for crosstalk-free quantum-gate operation in many-body nanofiber cavity QED systems*  
T. Keller, S. Kikura, R. Asaoka, Y. Suzuki, Y. Tokunaga, and T. Aoki 
[arXiv:2509.08408 (2025).](http://arxiv.org/abs/2509.08408)

Usage examples and the code for reproducing the paper data can be found in the notebook AddressingRequirements_arXiv.ipynb

## Usage
Specifying the total number of atoms N and the system type creates a cavity setup as detailed in the paper. Atoms targeted for performing a quantum gate are denoted by "x" and non-targeted atoms by "o" in the visualization. Photon frequency ωp is set on resonance with the bare cavity frequency ωc by default, which for simplicity is ωp = ωc = 0. All frequencies given in units of $2\pi$ x MHz.

    sys = System(N = 6, mode = :remote)
    Remote Gate | N = 6 | [xoo]--[xoo] | ωp = 0.0

By default we assume Cesium atoms with qubit splitting ωq, spontaneous emission rate 2γ and a transition coupling them to the cavity on resonance ωa = 0.

    atom = Atom()
    Atom | Target: false | ωq = 9192.6 | γ = 2.6 | g = 0.0 | ωa = 0.0

Create an empty cavity with default decay rate κr and mirror loss and transmission loss rates κm and κt respectively.

    cav = Cavity()
    Cavity [] | κr/κ = 0.9259 | (κr,κm,κt) = (2.50,0.10,0.10) | ωc = 0.0

Local and remote gate fidelities are calculated numerically either via the approximative superposition input...

    sys = System(N = 2, mode = :local)
    coupling = 7.8
    set_coupling(sys,coupling)
    input = normalize(ones(2^sys.N))
    target = CZ_unitary(input,sys)
    fidelity(local_gate(density_matrix(input),sys),target)

... or via the Choi operator

    sys = System(N = 2, mode = :remote)
    coupling = 7.8
    set_coupling(sys,coupling)
    input_choi = choi_state(sys.N)
    target_choi = choi(CZ_unitary(nothing,sys))*input_choi
    fidelity(remote_gate(nothing,sys,choi=true),target_choi,choi=true)

Similarly for gate success probabilities

    success_probability(sys)
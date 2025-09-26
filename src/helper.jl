const c = 2.99792458e8
const ħ = 1.05457182e-34
const kb = 1.380649e-23
const μ0 = 4π*1e-7
const ϵ0 = 8.8541878128e-12
const η0 = sqrt(μ0/ϵ0)
const a0 = 5.29177210544e-11
const e = 1.602176634e-19
const μe = 9.1093837139e-31
const amu = 1.66053906892e-27

function save_data(fname,labels,data)
    open(fname, "w") do file
        writedlm(file, labels)
        writedlm(file, data)
    end
end
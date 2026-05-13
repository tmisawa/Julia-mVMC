using MVMCExpertModeParsers: ExpertModeData

"""
Helpers for unit tests: electron-number vectors and simple moves.

This file intentionally stays lightweight and dependency-free so it can be
included from multiple unit test files without pulling in heavy machinery.
"""

"""
    make_ele_num(nsite; up_sites=[], down_sites=[])

Create an electron-occupation vector `ele_num` of length `2*nsite`.

Layout matches C/J conventions used in MVMCOptimizers:
- `ele_num[ri + 1]`             : up-spin occupation at site `ri` (0-based ri)
- `ele_num[nsite + ri + 1]`     : down-spin occupation at site `ri`
"""
function make_ele_num(
    nsite::Int;
    up_sites::Vector{Int} = Int[],
    down_sites::Vector{Int} = Int[],
)
    ele_num = zeros(Int, 2 * nsite)
    for ri in up_sites
        ele_num[ri + 1] += 1
    end
    for ri in down_sites
        ele_num[nsite + ri + 1] += 1
    end
    return ele_num
end

"""
    apply_hop(ele_num, ri, rj, s, nsite) -> new_ele_num

Apply a single hop move (ri -> rj) for spin `s` (0: up, 1: down) and return a
new electron-occupation vector.
"""
function apply_hop(ele_num::Vector{Int}, ri::Int, rj::Int, s::Int, nsite::Int)
    new_ele_num = copy(ele_num)
    new_ele_num[ri + s * nsite + 1] -= 1
    new_ele_num[rj + s * nsite + 1] += 1
    return new_ele_num
end

"""
    nsite(data)

Convenience accessor for tests that construct minimal `ExpertModeData`.
"""
nsite(data::ExpertModeData) = data.modpara.nsite


"""
Test suite for MVMCExpertModeParsers.jl
"""

using Test
using MVMCExpertModeParsers

# Include test files
include("test_parsers.jl")
include("test_green_two_ex_parser.jl")
include("test_integration.jl")
include("test_validation.jl")
include("test_utils.jl")
include("test_parse_expert_mode_files.jl")
include("test_parameter_initialization.jl")
include("test_read_input_parameters.jl")
include("test_read_input_parameters_rbm_layout.jl")
include("test_qp_weight.jl")
include("test_orbital_qptrans_utils.jl")
include("test_sfmt_compatibility.jl")
include("test_trans_parser_spin_indices.jl")
include("test_parameter_init_complexflag_rbm.jl")

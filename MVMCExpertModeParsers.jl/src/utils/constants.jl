"""
Constants for Expert Mode Parsing

Constants matching the C implementation and defining file formats.
"""

# File naming constants
const D_FILENAME_MAX = 256
const MAX_LINE_LENGTH = 1024

# Precision constants
const C_FLOAT_PRECISION = "%.18e"
const C_COMPLEX_FORMAT = "% .18e % .18e"
const C_INT_FORMAT = "%d"

# File extensions
const DEF_EXTENSION = ".def"
const NAMELIST_FILE = "namelist.def"

# Keyword definitions matching C implementation
const MVMC_KEYWORDS = Dict{String,String}(
    "ModPara" => "modpara.def",
    "LocSpin" => "locspn.def",
    "Trans" => "trans.def",
    "CoulombIntra" => "coulombintra.def",
    "CoulombInter" => "coulombinter.def",
    "Hund" => "hund.def",
    "Exchange" => "exchange.def",
    "PairHop" => "pairhop.def",
    "Gutzwiller" => "gutzwiller.def",
    "Jastrow" => "jastrow.def",
    "Orbital" => "orbital.def",
    "OrbitalParallel" => "orbitalparallel.def",
    "OrbitalAntiParallel" => "orbitalantiparallel.def",
    "OrbitalGeneral" => "orbitalgeneral.def",
    "QPTrans" => "qptrans.def",
    "OneBodyG" => "greenone.def",
    "TwoBodyG" => "greentwo.def",
    "TwoBodyGEx" => "greentwoex.def",
    "InterAll" => "interall.def",
    "OptTrans" => "opttrans.def",
    # Input parameter files (ReadInputParameters)
    "InGutzwiller" => "ingutzwiller.def",
    "InJastrow" => "injastrow.def",
    "InDH2" => "indh2.def",
    "InDH4" => "indh4.def",
    "InChargeRBM_PhysLayer" => "inchargerbmplyslayer.def",
    "InSpinRBM_PhysLayer" => "inspinrbmplyslayer.def",
    "InGeneralRBM_PhysLayer" => "ingeneralrbmplyslayer.def",
    "InChargeRBM_HiddenLayer" => "inchargerbmhiddenlayer.def",
    "InSpinRBM_HiddenLayer" => "inspinrbmhiddenlayer.def",
    "InGeneralRBM_HiddenLayer" => "ingeneralrbmhiddenlayer.def",
    "InChargeRBM_PhysHidden" => "inchargerbmphyshidden.def",
    "InSpinRBM_PhysHidden" => "inspinrbmphyshidden.def",
    "InGeneralRBM_PhysHidden" => "ingeneralrbmphyshidden.def",
    "InOrbital" => "inorbital.def",
    "InOrbitalAntiParallel" => "inorbitalantiparallel.def",
    "InOrbitalParallel" => "inorbitalparallel.def",
    "InOrbitalGeneral" => "inorbitalgeneral.def",
    "InOptTrans" => "inopttrans.def",
)

# Keyword indices matching C enum KWIdxInt
const KW_IDX_INT = Dict{String,Int}(
    "KWModPara" => 0,
    "KWLocSpin" => 1,
    "KWTrans" => 2,
    "KWCoulombIntra" => 3,
    "KWCoulombInter" => 4,
    "KWHund" => 5,
    "KWExchange" => 6,
    "KWPairHop" => 7,
    "KWGutzwiller" => 8,
    "KWJastrow" => 9,
    "KWDH2" => 10,
    "KWDH4" => 11,
    "KWOrbital" => 12,
    "KWOrbitalParallel" => 13,
    "KWOrbitalAntiParallel" => 14,
    "KWOrbitalGeneral" => 15,
    "KWQPTrans" => 16,
    "KWOneBodyG" => 17,
    "KWTwoBodyG" => 18,
    "KWTwoBodyGEx" => 19,
    "KWInterAll" => 20,
    "KWOptTrans" => 21,
    "KWTransSym" => 22,
    "KWBFRange" => 23,
    "KWBF" => 24,
    "KWIdxInt_end" => 25,
)

# Spin symbols
const SPIN_UP = :up
const SPIN_DOWN = :down
const SPIN_BOTH = :both

# Complex flags
const REAL_FLAG = 0
const COMPLEX_FLAG = 1

# Calculation modes
const PARAMETER_OPTIMIZATION = 0
const PHYSICS_CALCULATION = 1

# Lanczos modes
const LANCZOS_NONE = 0
const LANCZOS_ENERGY_ONLY = 1
const LANCZOS_GREEN_FUNCTIONS = 2

# Default values
const DEFAULT_NSITE = 0
const DEFAULT_NELEC = 0
const DEFAULT_NLOCSPIN = 0
const DEFAULT_VMC_CALC_MODE = 0
const DEFAULT_LANCZOS_MODE = 0
const DEFAULT_RND_SEED = 11272  # Expert Mode default (matches C implementation: readdef.c:1773)
const DEFAULT_NSR_OPT_ITR_STEP = 1000
const DEFAULT_NSR_OPT_ITR_SMP = 1000
const DEFAULT_NVMC_WARMUP = 1000
const DEFAULT_NVMC_INTERVAL = 1
const DEFAULT_NVMC_SAMPLE = 10000
const DEFAULT_DSR_OPT_RED_CUT = 1e-6
const DEFAULT_DSR_OPT_STA_DEL = 0.0
const DEFAULT_DSR_OPT_STEP_DT = 0.01
const DEFAULT_DSR_OPT_CG_TOL = 1e-6
const DEFAULT_NSR_OPT_CG_MAX_ITER = 1000
const DEFAULT_NSP_GAUSS_LEG = 1
const DEFAULT_NSP_STOT = 0
const DEFAULT_NMP_TRANS = 0
const DEFAULT_TWO_SZ = -1  # -1 means Sz is not conserved (FSZ mode), matching C default
const DEFAULT_N_DATA_IDX_START = 0
const DEFAULT_N_DATA_QTY_SMP = 1
const DEFAULT_C_DATA_FILE_HEAD = "zvo"
const DEFAULT_C_PARA_FILE_HEAD = "zqp"
const DEFAULT_N_FILE_FLUSH_INTERVAL = 1
const DEFAULT_COMPLEX_FLAG = 0
const DEFAULT_NNEURON = 0
const DEFAULT_NBLOCK_SIZE_RBM_RATIO = 8
const DEFAULT_N_ONE_BODY_G = 0
const DEFAULT_N_TWO_BODY_G = 0
const DEFAULT_N_TWO_BODY_G_EX = 0
const DEFAULT_NEX_UPDATE_PATH = 1
const DEFAULT_NSPLIT_SIZE = 1  # Default value for NSplitSize (from C SetDefaultValuesModPara)

# Error messages
const ERROR_FILE_NOT_FOUND = "File not found"
const ERROR_INVALID_FORMAT = "Invalid file format"
const ERROR_PARSE_ERROR = "Parse error"
const ERROR_VALIDATION_ERROR = "Validation error"
const ERROR_MISSING_PARAMETER = "Missing required parameter"
const ERROR_INVALID_VALUE = "Invalid parameter value"
const ERROR_INCONSISTENT_DATA = "Inconsistent data"

# Warning messages
const WARNING_DEFAULT_VALUE = "Using default value"
const WARNING_IGNORED_LINE = "Ignored line"
const WARNING_DEPRECATED_FORMAT = "Deprecated format"
const WARNING_APPROXIMATION = "Using approximation"

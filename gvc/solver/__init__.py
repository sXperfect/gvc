from .nn import NNSolver

SOLVERS = {
    "nn": NNSolver
}

AVAIL_SOLVERS = [k for k in SOLVERS.keys()]
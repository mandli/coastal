#!/usr/bin/env python
# encoding: utf-8

r"""Benchmark 1:  Solitary wave on a simple beach

Description:

"""

import sys

import numpy

import clawpack.pyclaw.util as util

def friction_source(solver, state, dt, manning=0.025, TOLERANCE=1e-30):
    r"""Compute the loss of momenta due to friction

    Solves:

        (hu)_t = -C_f |u|

    where C_f is determined via Manning's n formula.
    """
    g = state.problem_data['g']
    dry_tolerance = state.problem_data['dry_tolerance']
    
    if manning > TOLERANCE:
        wet_indices = numpy.nonzero(state.q[0,:] > dry_tolerance)
        gamma  = g * manning**2 / state.q[0, wet_indices]**(7.0 / 3.0)  \
                        * state.q[1, wet_indices]
        state.q[1, wet_indices] = state.q[1, wet_indices] / (1.0 + dt * gamma)


def solitary_wave(use_petsc=False, outdir="./_output"):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw


    # Physical/Algorithmic parameters
    g = 1.0
    sea_level = 0.0
    dry_tolerance = 1e-3

    # Parameters from the benchmark
    d = 1.0                                     # Constant depth away from beach
    H = d * 0.019                               # Height of wave
    gamma = numpy.sqrt(3.0 * H / (4.0 * d))     #
    L = numpy.arccosh(numpy.sqrt(20.0)) / gamma # Distance to soliton
    beta = 1.0 / numpy.arctan(19.85)            # Slope of beach (angle)
    x_0 = d / numpy.tan(beta)                   # Toe of beach
    x_1 = x_0 + L                               #

    # Setup solver
    solver = pyclaw.ClawSolver1D()
    solver.order = 2
    solver.cfl_max = 1.0
    solver.cfl_desired = 0.9
    solver.num_eqn = 2
    solver.num_waves = 2
    solver.limiters = pyclaw.limiters.tvd.MC
    solver.fwave = True
    solver.source_split = 0
    # Simplified Riemann solver, probably should switch to more complex one
    import rp_geoclaw
    solver.rp = rp_geoclaw
    #solver.step_source = \
    #    lambda solver, state, dt:friction_source(solver, state, dt, manning=0.0)

    # Boundary conditions
    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap

    # Setup domain
    x = pyclaw.Dimension(-10.0, 5.0, 700, name='x')
    domain = pyclaw.Domain(x)

    # Setup initial condition and bathymetry
    state = pyclaw.State(domain, 2, 1)

    # Bathymetry
    state.aux[0, :] = -d
    # state.aux[0, :] =   s(x.centers < x_0) * numpy.ones(x.centers.shape) * -d

                      # + (x.centers > x_0) * numpy.ones(x.centers.shape) \
                      # * (numpy.tan(beta) + d / (x.nodes[-1] - x_0)) * ((x.centers - x_0) - d)

    # Quiescent initial state
    state.q[0, :] = numpy.maximum(0.0, sea_level - state.aux[0, :])
    state.q[1, :] = 0.0

    # Now add solitary wave
    # eta = H / numpy.cosh(gamma * (x.centers - x_1) / d)**2
    # state.q[0, :] += (state.q[0, :] > dry_tolerance) * eta
    # state.q[1, :] = state.q[0, :] * -numpy.sqrt(g / d) * eta

    # Gaussian hump
    sigma = 0.2
    eta = 0.1 * numpy.exp(-x.centers**2/sigma**2)
    state.q[0, :] += (state.q[0, :] > dry_tolerance) * eta

    # Parameters for use in the Riemann solver
    state.problem_data['g'] = g
    state.problem_data['dry_tolerance'] = dry_tolerance
    state.problem_data['earth_radius'] = 0.0
    state.problem_data['deg2rad'] = 0.0
    state.problem_data['mcapa'] = 0

    # Controller
    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state, domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.keep_copy = False
    if use_petsc:
        claw.output_format = 'petsc'
    else:
        claw.output_format = 'ascii'
    claw.write_aux_init = True
    claw.output_style = 1
    claw.num_output_times = 9
    claw.tfinal = 1.0
    # claw.output_style = 3
    # claw.num_output_times = 10

    # Pass on some info to the plotting script
    with open('plot.data', 'w') as plot_data:
        plot_data.write("%s %s\n" % (x.lower, x.upper))
        plot_data.write("%s\n" % (claw.output_format))
        plot_data.write("%s\n" % dry_tolerance)

    return claw


if __name__ == "__main__":
    if len(sys.argv) > 1:
        pass

    # Compile Riemann solver if needed
    FFLAGS = "-O3 -funroll-loops -finline-functions"
    # util.compile_library(['rp1_shallow_fwave.f90'], 'rp_hll_fwave', FFLAGS=FFLAGS)
    util.compile_library(['rp1_geoclaw.f90', 'geoclaw_riemann_utils.f'], 'rp_geoclaw', FFLAGS=FFLAGS)

    from clawpack.pyclaw.util import run_app_from_main
    import setplot
    output = run_app_from_main(solitary_wave, setplot.setplot)


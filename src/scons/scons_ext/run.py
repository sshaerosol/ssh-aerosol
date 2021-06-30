# -*- coding: utf-8 -*-
# Copyright (C) 2016, ENPC - EDF R&D
#     Author(s): Sylvain Doré
#
# This file is part of the air quality modeling system Polyphemus.
#
# Polyphemus is developed in the INRIA - ENPC joint project-team CLIME and in
# the ENPC - EDF R&D joint laboratory CEREA.
#
# Polyphemus is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Polyphemus is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# For more information, visit the Polyphemus web site:
#      http://cerea.enpc.fr/polyphemus/

import os
from SCons.Script import *
from SCons.Script.SConscript import global_exports

# Checks the version of SCons.
EnsureSConsVersion(1, 0, 0)

# SCons is able to deal with duplicated environment.
SetOption('warn', 'no-duplicate-environment')

# BUILD_TARGETS has to be cleared before adding our own targets.
# Reasons why 'BUILD_TARGETS' is directly managed:
#  - Several targets can be deduced from the command line options
#     Example: scons mpi=* openmp=* my_algo
#       Despite that a single target was given, namely 'my_algo', this line
#       builds all the possible parallelization with automatic executable
#       renaming into "my_algo-mpi", "my_algo-omp" and "my_algo-mpi-omp".
#  - No need to specify '.exe' when building on Windows
#     (Yes, under Windows, '.exe' has to be added to the target name in
#      the scons command line)
while len(BUILD_TARGETS) > 0:
    BUILD_TARGETS.pop()

# The command line targets in full path.
command_line_target = [os.path.join(Dir('#').abspath, p)
                       for p in COMMAND_LINE_TARGETS][:]


def quiet(is_quiet=None):
    """Sets SCons quiet mode.

    @type quiet: bool
    @param quiet: 'True'/'False' if SCons should be quiet or not, 'None' to only
    get the current state.

    @return: The SCons quiet mode before this call.
    """
    former_display_mode = Main.progress_display.print_it
    if quiet is not None:
        Main.progress_display.set_mode(not is_quiet)
    return not former_display_mode


def sconstruct(script, _cache=dict()):
    """Executes the SConstruct file 'script'."""

    if os.path.isfile(script + "/SConstruct"):
        script += "/SConstruct"
    script_path = os.path.abspath(script)
    if script_path in _cache:
        return _cache[script_path]

    global global_exports
    initial_exports = global_exports.copy()
    is_quiet = quiet()
    try:
        Export(src_dir = os.path.dirname(os.path.abspath(script)))
        if not is_quiet:
            print("[SCONSTRUCT]", script)
        targets = SConscript(script)
        if targets is None:
            targets = []
        if type(targets) is not list:
            targets = [targets]
        _cache[script_path] = targets
    finally:
        # It global exports have been made, then it should not be visible
        # from other SConstruct files:
        global_exports.clear()
        global_exports.update(initial_exports)
        quiet(is_quiet)

    return targets


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

import glob
import os
from SCons.Script import *

from . import run


class Spack:
    def __init__(self, utils):
        self.utils = utils

    def dependency(self, path, exclude_dependency, build_dir=None):
        path = os.path.abspath(path)
        species_file = glob.glob(os.path.join(path, "*.species")) + \
                       glob.glob(os.path.join(path, "species"))
        reactions_file = glob.glob(os.path.join(path, "*.reactions")) + \
                         glob.glob(os.path.join(path, "reactions"))
        spack_config = glob.glob(os.path.join(path, "spack_config"))
        if spack_config:
            spack_config = spack_config[0]

        if not species_file or not reactions_file:
            if spack_config:
                raise Exception, "In \"" + path + "\":\n"\
                                 "Spack needs a '.species' file "\
                                 "and a '.reactions' file to operate"
            return []
        if len(species_file) > 1:
            if spack_config:
                raise Exception, "In \"" + path + "\":\n"\
                                 "Several species file given to Spack, "\
                                 "only one is expected"
            return []
        if len(reactions_file) > 1:
            if spack_config:
                raise Exception, "In \"" + path + "\":\n"\
                                 "Several reactions file given to Spack, "\
                                 "only one is expected"
            return []
        if not spack_config:
                raise Exception, "In \"" + path + "\":\n"\
                        "There is a species file and reactions file, "\
                        "but Spack needs a 'spack_config' file to operate"


        # Creates a vanilla environment since Spack does not depends on
        # environment specifics and can be run in the source tree instead
        # of the build directory.
        spack_env = self.utils.create_env()

        spack_fortran_output = []

        for filename in ["dimensions.f90", "dratedc.f90", "fexchem.f90",
                         "jacdchemdc.f90", "kinetic.f90", "rates.f90",
                         "LU_decompose.f90", "LU_solve.f90", "LU_solve_tr.f90"]:
                         # use in two-step time numerical solver
                         # "fexloss.f90", "fexprod.f90"]:

            fortran_file = os.path.join(path, filename)
            spack_fortran_output.append(fortran_file)

        spack_config_output = os.path.join(path, "species.spack.dat")

        spack_exe = os.path.join(self.utils.polyphemus_path,
                                "include/spack/src/generate_chem_files.sh")
        spack_arg = [os.path.basename(p)
                     for p in species_file + reactions_file]
        spack_input = [spack_config] + species_file + reactions_file

        # Generates the files with Spack.
        spack_target = spack_env.Command(spack_fortran_output +
                                         [spack_config_output],
                                         spack_input,
                                         [ "cd " + path + " && " +
                                          " ".join([spack_exe] + spack_arg)])

        # Automatically builds Spack if needed.
        spack_dir = os.path.dirname(spack_exe)
        run.sconstruct(spack_dir + "/SConstruct")
        Depends(spack_target, Dir(spack_dir))

        # Copies the generated files into the current build directory.
        copy_fortran = []
        for f in spack_fortran_output:
            copy_fortran.append(Copy(self.utils.rebase_dir(build_dir, f), f))
        fortran_target = spack_env.Command(self.utils.rebase_dir(build_dir,
                                                                 spack_fortran_output),
                                           spack_target,
                                           copy_fortran)

        # Avoids building Spack files from the source tree, we have copies
        # in the build directories.
        exclude_dependency += spack_fortran_output
#        return [fortran_target]
        return fortran_target # YK

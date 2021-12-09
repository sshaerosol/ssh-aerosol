# -*- coding: utf-8 -*-
# Copyright (C) 2007-2016, ENPC - INRIA - EDF R&D
#     Author(s): Sylvain Doré, Vivien Mallet
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

import subprocess, glob, os, re, sys
import distutils.sysconfig
from SCons.Errors import UserError
from SCons.Script import *
from SCons.Script.SConscript import global_exports

from . import FlagDict, Spack, load_profile, list_profile
from . import getstatusoutput


class Utils:
    """Helps to configure SCons based on the user defined build variables."""

    # Global scons extensions accessible through 'Utils' instances.
    from . import command_line_target
    from . import run


    def __init__(self, user_variables, exported_variables):
        """Constructs an instance of 'Utils'.

        @type user_variables: dictionary
        @param user_variables: The user variables used for the build.

        @type _exported_variables: dictionary
        @param _exported_variables: The global variables used for the build,
        globally exported with the SCons 'Export' function.
        """
        self._user_variables = user_variables
        self._exported_variables = exported_variables
        self.profile = FlagDict()

        #=== Shared variables
        self.src_dir = self.create_variable("src_dir", Dir('#').abspath)

        polyphemus_path = self.create_variable("polyphemus_path", None)
        polyphemus_path = ARGUMENTS.get("polyphemus", polyphemus_path)
        if polyphemus_path is None:
            polyphemus_path = Dir('#').abspath
        elif not os.path.isabs(polyphemus_path):
            polyphemus_path = os.path.join(self.src_dir, polyphemus_path)
        polyphemus_path = os.path.abspath(polyphemus_path)

        if not os.path.isdir(polyphemus_path):
            raise Exception("The Polyphemus path \"" + polyphemus_path \
                            + "\" does not appear to be a valid path.")
        self.polyphemus_path = polyphemus_path

        self.scons_dir = self.create_variable("scons_dir", self.polyphemus_path)

        # Comment when spack is turned off.
        #=== Custom scons extensions.
        self.spack = Spack(self)


    def create_variable(self, name, default):
        """Returns the value of the variable whose name is in 'name'.

        The variable 'name' is looked up in the user defined variables, then in
        the variables exported with the SCons 'Export' function. If the variable
        is not found or is 'None', then 'default' is returned.
        """
        value = None
        if self._user_variables is not None \
            and name in self._user_variables:
            value = self._user_variables[name]
        elif self._exported_variables is not None \
            and name in self._exported_variables:
            value = self._exported_variables[name]
        if value is not None:
            return value
        else:
            return default


    def create_list(self, name, possible_values = None, default = None):
        """Returns a list from the variable 'name'.

        The list is constructed from the variable 'name'. Then the value
        of 'name' in the profile and in the command line arguments are appended
        in this order. If the list is empty, 'default' is used, otherwise the
        empty list is returned.

        If 'possible_values' is not None:
          - values are checked against it.
          - if '*' is found, then 'possible_values' is returned.
        """
        l = self.create_variable(name, [])
        if not SCons.Util.is_List(l):
            l = str(l)
        l = Split(l)
        l += Split(self.profile.get(name, ""))
        l += Split(ARGUMENTS.get(name, ""))
        if not l:
            if default is None:
                default = []
            elif type(default) in (str, unicode):
                default = [default]
            l = default

        if possible_values:
            if '*' in l:
                return possible_values
            for v in l:
                self.assert_valid_argument(name, v, possible_values)

        return l


    def create_path_list(self, name, default = None):
        """Returns a path list from the variable 'name'.

        Relative path are interpreted from the current directory, and, if they
        don't exist from the current directory, from the Polyphemus path.
        """
        path_list = []
        for path in self.create_list(name):
            if os.path.isdir(path):
                path_list.append(path)
            elif os.path.isdir(os.path.join(self.polyphemus_path, path)):
                path = os.path.join(self.polyphemus_path, path)
                path_list.append(path)
            else:
                raise Exception("Unable to find the '" + name + "' directory "\
                                "\"" + path + "\" (even in Polyphemus directory, \"" \
                                + self.polyphemus_path + "\").")
        return path_list


    def create_flag_string(self, name):
        """Returns a space separated string for the flag variable 'name'.

        The space separated string is constructed from the value returned by
        'create_list(name)'.
        """
        return " ".join(self.create_list(name))


    def assert_valid_argument(self, name, value, possible_value):
        if value not in possible_value:
            raise Exception("Unsupported option \"" + value \
                            + "\" for argument \"" + name \
                            + "\". Available options are: " \
                            + ", ".join(["\"" + x + "\"" for x in possible_value]) + ".")


    def add_argument(self, name, value_list = None):
        """Checks a command line argument value, sets a default value if needed.

        Checks whether the command line argument 'name' has been provided, and
        sets the argument to its default value if the user has not given any
        value.
        """
        if value_list is None:
            if not ARGUMENTS.has_key(name):
                raise Exception("The command line argument \"" + name + \
                      + "\" is required, but it was not provided.")
        else:
            ARGUMENTS[name] = ARGUMENTS.get(name, value_list[0])
            self.assert_valid_argument(name, ARGUMENTS[name], value_list)


    def debug_flag(self, name):
        """Returns the debug flags from the command line argument 'name'."""
        if ARGUMENTS[name] == "-1":
            return ""
        elif ARGUMENTS[name] == "0":
            # '-g' has no impact on performances.
            return "-O2 -g"
        elif ARGUMENTS[name] == "1":
            return "-g"
        elif ARGUMENTS[name] == "2":
            return "-O2 -g"


    def rebase_dir(self, base_dir, path_list):
        """Returns path starting from 'base_dir' instead of 'polyphemus_path'."""
        if base_dir is None:
            return path_list
        root_len = len(self.polyphemus_path)

        def translate_path(p):
 #           if type(p) == str or type(p) == unicode:
            if type(p) == str:
                 return base_dir + os.path.abspath(p)[root_len:]
            else:
                return p

        if type(path_list) == list:
            return list([translate_path(p) for p in path_list])
        else:
            return translate_path(path_list)


    def create_env(self):
        """Creates a generic Polyphemus SCons environment."""
        # F90PATH='.' is needed to have a correct build order for F90 files.
        env = Environment(ENV = os.environ, F90PATH='.')
        env.Replace(CONFIGURELOG = "#/.scons.log")
        SCons.Script.SConsignFile(self.scons_dir + '/.sconsign.dblite')
        return env


    def create_programs(self, generate_program = True):
        """Creates the SCons program targets."""

        #############
        # ARGUMENTS #
        #############

        self.add_argument("debug", ["-1", "0", "1", "2"])
        self.add_argument("debug_cpp", [ARGUMENTS["debug"], "-1", "0", "1", "2"])
        self.add_argument("debug_fortran",
                          [ARGUMENTS["debug"], "-1", "0", "1", "2"])
        self.add_argument("line", ["no", "yes"])
        self.add_argument("mode_cpp", ["strict", "permissive"])
        self.add_argument("mode_fortran", ["permissive", "strict"])
        self.add_argument("_", ["1", "0", "2"])
        self.add_argument("openmp", ["no", "yes"])
        self.add_argument("mpi", ["no", "yes"])
        self.add_argument("intel", ["no", "yes"])

        has_deprecated_debug = ARGUMENTS["debug"] != "-1" \
                            or ARGUMENTS["debug_cpp"] != "-1" \
                            or ARGUMENTS["debug_fortran"] != "-1"

        ########################
        # ENVIRONMENT CREATION #
        ########################

        env = self.create_env()

        if "LD_LIBRARY_PATH" in os.environ:
            env.Append(LIBPATH = os.environ["LD_LIBRARY_PATH"].split(os.pathsep))
        if "LIBRARY_PATH" in os.environ:
            env.Append(LIBPATH = os.environ["LIBRARY_PATH"].split(os.pathsep))
        if "CPATH" in os.environ:
            env.Append(CPPPATH = os.environ["CPATH"].split(os.pathsep))
        if "CPLUS_INCLUDE_PATH" in os.environ:
            env.Append(CPPPATH = os.environ["CPLUS_INCLUDE_PATH"].split(os.pathsep))

        if ARGUMENTS["line"] == "no":
            env.Replace(CCCOMSTR = "[C] $SOURCE")
            env.Replace(CXXCOMSTR = "[C++] $SOURCE")
            env.Replace(F77COMSTR = "[F77] $SOURCE")
            env.Replace(F77PPCOMSTR = "[F77-PP] $SOURCE")
            env.Replace(F90COMSTR = "[F90] $SOURCE")
            env.Replace(FORTRANCOMSTR = "[FORTRAN] $SOURCE")
            env.Replace(FORTRANPPCOMSTR = "[FORTRAN-PP] $SOURCE")
            env.Replace(LINKCOMSTR = "[Linking] $TARGET")


        ########################
        # COMPILERS AND LINKER #
        ########################

        c_compiler = self.create_variable("c_compiler", None)
        cpp_compiler = self.create_variable("cpp_compiler", None)
        fortran_compiler = self.create_variable("fortran_compiler", None)
        linker = self.create_variable("linker", "$CXX")

        if ARGUMENTS["intel"] == "yes":
            c_compiler = "icc"
            cpp_compiler = "icpc"
            fortran_compiler = "ifort"


        if ARGUMENTS["mpi"] == "yes":
            if ARGUMENTS["intel"] == "yes":
                cpp_compiler = "mpic++"
                c_compiler = "mpicc"
                fortran_compiler = "mpif90"
            else:
                for cpp_compiler in ['mpiCC', 'mpicxx', 'mpic++']:
                    if WhereIs(cpp_compiler) != None:
                        break
                    else:
                        raise Exception("Unable to find a MPI compiler.")
                    c_compiler = "mpicc"
                    fortran_compiler = "mpif90"
                    linker = cpp_compiler

        # The compilers and the linker can be overridden from the command line.
        c_compiler = ARGUMENTS.get("c", c_compiler)
        cpp_compiler = ARGUMENTS.get("cpp", cpp_compiler)
        fortran_compiler = ARGUMENTS.get("fortran", fortran_compiler)
        linker = ARGUMENTS.get("link", linker)

        fortran_env = ['F77', 'F90', 'F95', 'FORTRAN']
        if c_compiler is not None:
            env.Replace(CC = c_compiler)
        if cpp_compiler is not None:
            env.Replace(CXX = cpp_compiler)
        if fortran_compiler is not None:
            for v in fortran_env:
                env.Replace(**{v : fortran_compiler})
        if linker is not None:
            env.Replace(LINK = linker)

        if 'toolchain_suffix' in ARGUMENTS:
            suffix = ARGUMENTS['toolchain_suffix']
            suffix_env_name = "CUSTOM_TOOLCHAIN_SUFFIX"
            env[suffix_env_name] = suffix
            for tool in ['CC', 'CXX', 'LINK'] + fortran_env:
                current = env.subst('$' + tool)
                if current and not current.endswith(suffix):
                    env.Append(**{tool : "$" + suffix_env_name})

        ################
        # PROFILE LOAD #
        ################

        self.profile = FlagDict()
        if not env.GetOption('clean'):
            if has_deprecated_debug:
                print("""
[WARNING] === PLEASE READ: ===
[WARNING] 'debug=' is deprecated, please use 'profile=' instead. It gives access
to highly recommended debugging and optimization options.
[WARNING] ====================

""")
            elif 'profile' not in ARGUMENTS:
                profile_compiler, profile_dict = list_profile(env)
                if not profile_dict:
                    print("[WARNING] No profile available for compiler '{0}'.")\
                        .format(profile_compiler)
                else:
#                     raise UserError, """
# [ERROR] No profile given, please use 'profile=NAME'...
# [ERROR] (You can use 'profile=?' to get the list of available profile)"""
                    # Set the default profile to fast (YK: 2018/03/14)
                    ARGUMENTS["profile"] = 'fast'
                    load_profile(self, env, ARGUMENTS["profile"])
            else:
                load_profile(self, env, ARGUMENTS["profile"])


        ####################
        # COMPILER OPTIONS #
        ####################

        #=== Preprocessor options.
        preprocessor_defines = []
        preprocessor_defines_str = self.create_list("preprocessor_defines")
        for define in preprocessor_defines_str:
            preprocessor_defines.append(define.split("=", 1))
        if ARGUMENTS["_"] == "1":
            preprocessor_defines.append("POLYPHEMUS_SINGLE_UNDERSCORE")
        elif ARGUMENTS["_"] == "2":
            preprocessor_defines.append("POLYPHEMUS_DOUBLE_UNDERSCORE")

        #=== C++-specific compilation options.
        cpp_compilation_option = self.debug_flag("debug_cpp")
        # In case of GNU compilers, a few options are added.
        if "gcc" in env["CC"] and ARGUMENTS["mode_cpp"] == "strict":
            cpp_compilation_option += " -Wall -ansi -pedantic -Wno-unused" \
                + " -Wno-parentheses"
 
        #=== Fortran-specific compilation options.
        fortran_compilation_option = self.debug_flag("debug_fortran")
        if "gfortran" in env["FORTRAN"]:
            if ARGUMENTS["mode_fortran"] == "strict":
                fortran_compilation_option += " -Wall -pedantic"
            # Adds preprocessor to Fortran compilation.
            fortran_compilation_option += " -cpp"

        #=== OpenMP options.
        if ARGUMENTS["openmp"] == "yes":
            if "g++" in env["CXX"]:
                cpp_compilation_option += " -fopenmp"
                preprocessor_defines.append("POLYPHEMUS_PARALLEL_WITH_OPENMP")
            elif "icpc" in env["CXX"]:
                cpp_compilation_option += " -openmp"
                preprocessor_defines.append("POLYPHEMUS_PARALLEL_WITH_OPENMP")
            elif ARGUMENTS.has_key("flag_openmp"):
                cpp_compilation_option += " " + ARGUMENTS["flag_openmp"]
                preprocessor_defines.append("POLYPHEMUS_PARALLEL_WITH_OPENMP")
            else:
                print("[WARNING]: No openMP parallelization. Please use the" \
                    " option 'flag_openmp' to indicate the openMP compiling" \
                    " option to your C++ compiler.")

            preprocessor_defines.append("BZ_THREADSAFE")

            if "gfortran" in env["FORTRAN"]:
                fortran_compilation_option += " -fopenmp"
            elif "ifort" in env["FORTRAN"]:
                fortran_compilation_option += " -openmp"
            elif ARGUMENTS.has_key("flag_openmp"):
                fortran_compilation_option += " " + ARGUMENTS["flag_openmp"]
            else:
                print("[WARNING]: No openMP parallelization. Please use the" \
                    " option flag_openmp to indicate the openMP compiling" \
                    " option suitable to your FORTRAN compiler.")

        #=== MPI options.
        if ARGUMENTS["mpi"] == "yes":
            preprocessor_defines.append("POLYPHEMUS_PARALLEL_WITH_MPI")

        #=== FastJX options.
        # For enabling FastJX when it has been installed.
        fastjx_file = os.path.join(self.polyphemus_path, "include/fastJX/fastJX.f")
        if os.path.isfile(fastjx_file):
            preprocessor_defines.append("POLYPHEMUS_FASTJX")

        #=== User provided options.
        # Most compilers will give them precedence since they are appended last.
        flag_compiler = self.create_flag_string("flag_compiler")
        if flag_compiler:
            cpp_compilation_option += " " + flag_compiler
            fortran_compilation_option += " " + flag_compiler

        flag_cpp = self.create_flag_string("flag_cpp")
        if flag_cpp:
            cpp_compilation_option += " " + flag_cpp

        flag_fortran = self.create_flag_string("flag_fortran")
        if flag_fortran != "":
            fortran_compilation_option += " " + flag_fortran

        env.Append(CPPDEFINES=preprocessor_defines)
        env.Replace(CCFLAGS = cpp_compilation_option)
        env.Replace(F77FLAGS = fortran_compilation_option)
        env.Replace(FORTRANFLAGS = fortran_compilation_option)
        env.Replace(F90FLAGS = fortran_compilation_option) # YK

        ################
        # INCLUDE PATH #
        ################

        #=== User-defined paths for includes.
        # (It is also used for finding sources and SConstruct dependencies.)
        include_path = self.create_path_list("include_path")

        #=== Include search path for the compiler.
        include_search_path = self.create_path_list("include_search_path")
        include_search_path += include_path

        library_list = self.create_list("library_list")
        library_list += Split(ARGUMENTS.get("library_list", ""))
        if "python" in library_list:
            include_search_path.append(distutils.sysconfig.get_python_inc())

        env.Append(CPPPATH = include_search_path)
        env.Append(F77PATH = include_search_path)
        env.Append(FORTRANPATH = include_search_path)
        include_search_path.append("/usr/include/") # YK
        env.Append(F90PATH = include_search_path) # YK

        #############
        # LIBRARIES #
        #############

        #=== Sets the library search path list.
        library_path = self.create_path_list("library_path")
        env.Append(LIBPATH = library_path)

        #=== Adds "common" libraries.
        # TODO: only really used libraries should be linked.
        for library in ["ifcore", "imf", "svml", "intlc", "dl",
                        "m", "stdc++", "blas", "atlas", "lapack", "g2c",
                        "gfortran", "gslcblas", "blitz", "netcdff", "netcdf",
                        "netcdf_c++"]:
            if library not in library_list:
                library_list += [library]

        #=== Checks which libraries are available.
        conf = Configure(env, conf_dir=self.scons_dir + "/.sconf_temp")
        for library in library_list:
            if library == "python":
                python_version = distutils.sysconfig.get_python_version()
                library = "python" + python_version
            conf.CheckLib(library)
        env = conf.Finish()


        ##################
        # LINKER OPTIONS #
        ##################

        link_options = ""

        #=== OpenMP options.
        if ARGUMENTS["openmp"] == "yes":
            if "g++" in env["CXX"]:
                link_options += " -fopenmp"
            elif "icpc" in env["CXX"]:
                link_options += " -openmp"
            elif ARGUMENTS.has_key("flag_openmp"):
                link_options += " " + ARGUMENTS["flag_openmp"]
            else:
                print("[WARNING]: No openMP parallelization. Please, use the " \
                      "options flag_openmp to add the appropriate openMP " \
                      "linking option.")

        #=== User provided options.
        # Most linkers will give them precedence since they are appended last.
        flag_link = self.create_flag_string("flag_link")
        if flag_link != "":
            link_options += " " + flag_link

        env.Replace(LINKFLAGS = link_options)


        ########################
        # BUILD OUTPUT OPTIONS #
        ########################

        #=== Build output directory.
        build_dir = self.create_variable("build_dir",
                                         self.polyphemus_path + "/.build/")

        build_flavor = self.create_list("build_flavor")

        # Spack is not complied because of this flavor added (YK).
        # if ARGUMENTS["mpi"] == "yes":
        #     build_flavor.append("mpi")
        # if ARGUMENTS["openmp"] == "yes":
        #     build_flavor.append("omp")

        build_flavor_str = '-'.join(build_flavor)

        if build_flavor_str:
            build_dir += build_flavor_str
        else:
            build_dir += "default"

        # Builds everything inside build_dir:
        env.VariantDir(build_dir, self.polyphemus_path, duplicate=0)

        #=== Executable suffix.
        program_suffix_str = ""
        # Executable suffix contains the build flavor:
        if build_flavor_str:
            program_suffix_str += '-' + build_flavor_str
        program_suffix = self.create_list("program_suffix")
        if program_suffix:
            program_suffix_str += '-' + '-'.join(program_suffix)
        if program_suffix_str:
            env["PROGSUFFIX"] = program_suffix_str + env["PROGSUFFIX"]


        ############
        # PROGRAMS #
        ############

        #=== The targets to be built.
        target_list = self.create_list("target_list")
        if not target_list:
            target_list = glob.glob(self.src_dir + "/*.cpp")
            target_list += glob.glob(self.src_dir + "/ssh-aerosol.f90") # YK

        # In case there is a list of targets to be excluded.
        exclude_target = self.create_list("exclude_target")
        for target in target_list[:]:
            if target in exclude_target or target[:-4] in exclude_target:
                target_list.remove(target)

        #=== Directory dependencies.
        dir_dependencies = []
        for path in include_path:
            if os.path.isfile(path + "/SConstruct"):
                if os.path.abspath(path) !=  os.path.abspath(self.src_dir):
                    dir_dependencies += self.run.sconstruct(path)

        #=== Source dependencies.
        src_dependencies = self.create_list("src_dependencies")
        exclude_dependency = self.create_list("exclude_dependency")

        for path in include_path:
            src_dependencies += glob.glob(os.path.join(path, "*.[fFcC]"))
            src_dependencies += glob.glob(os.path.join(path, "*.f90"))
            src_dependencies += glob.glob(os.path.join(path, "*.F90"))
            # Comment when spack is turned off.
            src_dependencies += self.spack.dependency(path, exclude_dependency,
                                                      build_dir) # YK

        # In case there is a list of dependencies to be excluded.
        filtered_dependencies = []
        exclude_dependency = []
        # exclude_dependency = [ re.compile(p) for p in exclude_dependency ][:]
        # regex = re.compile(r'.build')
        for dependency in src_dependencies:
            regex = re.compile(str(dependency))
            if not exclude_dependency:
                filtered_dependencies.append(dependency)
            else:
                is_excluded = False
                str_dep = str(dependency)
                for expression in exclude_dependency:
                    if regex.search(expression) is not None:
                        is_excluded = True
                        break
                if str_dep in target_list:
                    is_excluded = True
                if not is_excluded:
                    filtered_dependencies.append(dependency)

        # for dependency in src_dependencies:
        #     str_dep = str(dependency)
        #     for expression in exclude_dependency:
        #         print "expression", expression, str_dep
        #         if expression.search(str_dep) is not None:
        #             break
        #     else:
        #         if str_dep in target_list:
        #             break
        #         else:
        #             print "*** appended dependency", dependency
        #             filtered_dependencies.append(dependency) # YK

            
        src_dependencies = filtered_dependencies

        #=== Programs / shared library creation.
        # Rebasing dependencies and target into 'build_dir':
        # (to call BuildVariant is not enough, this is a SCons strangeness)
        src_dependencies = self.rebase_dir(build_dir, src_dependencies)

        # Informs SCons on the targets to be built. It is assumed that all
        # targets have the same dependencies.
        for target in target_list:
            program_name = '.'.join(os.path.splitext(target)[:-1])
            program_dependencies = [self.rebase_dir(build_dir, target)] \
                                    + src_dependencies + dir_dependencies
            # Either we generate the program
            if generate_program:

                if (('ssh-aerosol.f90' in target) and
                    (ARGUMENTS["intel"] == "yes")):
                    env.Replace(LINK = fortran_compiler) # YK
                
                program = env.Program(program_name, program_dependencies)
                if program_name in self.command_line_target:
                    BUILD_TARGETS.append(program_name + env["PROGSUFFIX"])

                # In case another SConstruct wants to depend on this source
                # directory:
                # (It is used to build Spack as needed)
                Depends(Dir(self.src_dir), program)

            # Or we generate the shared library for ssh-aerosol
            else:
                if 'ssh-aerosol.f90' in target:
                    env.SharedLibrary('ssh-aerosol.so', program_dependencies)

        return env

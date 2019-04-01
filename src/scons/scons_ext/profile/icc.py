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

from scons_ext import FlagDict
from SCons.Script import ARGUMENTS

# There are changes at each GCC version, read with caution.
# For advanced debugging functionalities, see:
# https://gcc.gnu.org/onlinedocs/gcc/Instrumentation-Options.html

# NOTE for the sanitizer:
# The env variables ASAN_OPTIONS, LSAN_OPTIONS and TSAN_OPTIONS can define
# some of the sanitizer behavior.


def _warning(utils):
    p = utils.profile

    # p.flag_cpp += " -Wall -Wextra -Wundef"
    # -Wno-unused -Wno-parentheses are added to avoid too many warnings from gcc (YK).
    p.flag_cpp += " -Wall"

    # GFortran gets '-Wall' minus '-Wconversion -Wno-tabs -Wunused'
    # (with '-Wall' GFortran would become very chatty...)
    # p.flag_fortran += " -Waliasing -Wampersand -Wsurprising -Wc-binding-type" \
    #     " -Wintrinsics-std  -Wintrinsic-shadow -Wline-truncation" \
    #     " -Wtarget-lifetime -Wreal-q-constant"
    # -Wc-binding-type is not available in old gfortran (YK).
    p.flag_fortran += "  "

    # Colors the compiler output.
    # ('auto' restricts coloring to terminal output.)
    # -fdiagnostics-color=auto is not available in old version (YK).
    # p.flag_compiler += " -fdiagnostics-color=auto"

    return p


def fast_portable(utils):
    """Optimization restricted to be compatible with older CPUs.

    If compatibility with older computers is not needed, use 'fast' instead.
    """
    p = _warning(utils)
    p.flag_compiler += " -O2 -g3"
    p.flag_compiler += " -fPIC"
    return p


def fast(utils):
    """Optimization using full capacities of this computer's CPU.

    The resulting executable will take advantage of this computer's CPU and
    won't work on computers with older CPUs. Newer CPUs should be backward
    compatible.
    """
    p = fast_portable(utils)
    return p


def debug(utils):
    """Debugging options.

    This is the configuration to use when developing. If it gets too slow,
    there are two solutions:
     - add a command line argument to enable/disable slower functionalities,
     - migrate slower functionalities to another profile.

    Command line arguments to alter this profile:
     undefined_behavior=(abort|log)  default: abort
         When an undefined behavior is encountered (includes invalid numerical
         operations among miscellaneous checks), abort the process
         or just display a log message
     float_div_by_zero=(abort|log|ignore)  default: abort
         When a float division by zero is encountered, abort the process,
         just display a log message or simply ignore.
     float_overflow=(abort|ignore)  default: abort
         When a float overflow is encountered, abort the process,
         or simply ignore.
    optim=(none|light|fast)  default: light
        Use 'none' to disable any optimization, 'light' to enable debugger
        friendly optimization (the default) and 'fast' for maximum optimization.
    """
    utils.add_argument("undefined_behavior", ["abort", "log"])
    utils.add_argument("float_div_by_zero",  ["abort", "log", "ignore"])
    utils.add_argument("float_overflow",     ["abort", "ignore"])
    utils.add_argument("optim",              ["none", "light", "fast"])

    ############
    # Warnings #
    ############

    p = _warning(utils)
    # An additional _warning that gives performance hints.
    #p.flag_cpp += " -Weffc++"

    ##########
    # Macros #
    ##########

    p.preprocessor_defines += " " + " ".join([
        # Activates glibc buffer overflow checking, notably for:
        # memcpy, mempcpy, memmove, memset, stpcpy, strcpy, strncpy,
        # strcat, strncat, sprintf, snprintf, vsprintf, vsnprintf,
        # and gets.
        "_FORTIFY_SOURCE=2",

        # libstdc++ debug mode.
        "_GLIBCXX_DEBUG",
        # Enforces undefined/forbidden behaviors from the standard.
        "_GLIBCXX_DEBUG_PEDANTIC",

        # Blitz debug mode.
        "BZ_DEBUG",

        # Polyphemus related macros.
        "TALOS_DEBUG",
        "SELDONDATA_DEBUG_LEVEL_4",
        "SELDON_DEBUG_LEVEL_4",
        "VERDANDI_DEBUG_LEVEL_4",
        "SELDON_WITH_ABORT",
        "VERDANDI_WITH_ABORT",
        "OPS_WITH_ABORT"
    ])

    ######################
    # GCC built-in debug #
    ######################
    if ARGUMENTS["optim"] == "none":
        p.flag_compiler +=  " -O0"
    elif ARGUMENTS["optim"] == "light":
        # -Og is optimizing while prioritizing debugging
        #     (sets -fno-omit-frame-pointer -fno-inline etc.)
        p.flag_compiler +=  " -Og"
    elif ARGUMENTS["optim"] == "fast":
        # Using portable optimization since it will be slow anyway.
        fast_portable(utils)

    p.flag_compiler += " -march=native"

    p.flag_compiler +=  " " + " ".join([
        # -g3 is -g2 (the default) with macro symbols added.
        "-g3",

        # Adds stack memory guards to detect buffer overflows.
        "-fstack-protector-all",

        # Undefined behavior sanitizer, aka ubsan.
        # (This flag enables common ubsan functionalities.)
        # It is not working in old version (YK).
       # "-fsanitize=undefined",

        # Floating-point to integer conversion checking.
        # "-fsanitize=float-cast-overflow",
            # This option enables floating-point type to integer conversion checking.
            # We check that the result of the conversion does not overflow.  Unlike
            # other similar options, -fsanitize=float-cast-overflow is not enabled by
            # -fsanitize=undefined because this option does not work well with the
            # useful invalid floating point exception, aka "FE_INVALID".
            # The macro TALOS_DEBUG enables "FE_INVALID", please define the macro
            # 'TALOS_NO_FE_INVALID' if using '-fsanitize=float-cast-overflow'.
    ])

    if ARGUMENTS["float_overflow"] != "abort":
        p.preprocessor_defines += " TALOS_NO_FE_OVERFLOW"

    abort_list = []
    if ARGUMENTS["undefined_behavior"] == "abort":
        abort_list.append("undefined")

    # This is not available in old version (YK).
    # if ARGUMENTS["float_div_by_zero"] == "log":
    #     p.flag_compiler += " -fsanitize=float-divide-by-zero"
    if ARGUMENTS["float_div_by_zero"] != "abort":
        # 'TALOS_NO_FE_DIVBYZERO' prevents Talos to abort the program on the first
        # division by zero. We prefer to rely on the  sanitizer to do this.
        p.preprocessor_defines += " TALOS_NO_FE_DIVBYZERO"

    # TO FIX: SANITIZER limitation at least in gcc-5.3
    #         cannot produce core dump on float division by zero
    #         (also, cannot find an handle to break on with gdb)
#         abort_list.append("float-divide-by-zero")
# -fno-sanitize-recover= is not working (YK).
#    if abort_list:
#        p.flag_compiler += " -fno-sanitize-recover=" + ",".join(abort_list)


    ###########################
    # GFortran built-in debug #
    ###########################
    p.flag_fortran += " " + " ".join([
        # All sort of runtime check - includes bounds checking.

        # The run-time error 'recursive call to nonrecursive procedure' appears
        #        "-fcheck=all", (YK)
        "-fcheck=array-temps,bounds,do,mem,pointer",
        # Equivalent to "implicit none" in every procedure.
        #-fimplicit-none
        # Display the stack on fatal error
        "-fbacktrace",
    ])

    # GFortran can trap invalid operations and float overflow. It can also trap
    # float division by zero, but we rely on the sanitizer for this.
    if ARGUMENTS["float_overflow"] == "abort":
        p.flag_fortran += " -ffpe-trap=invalid,overflow"
    else:
        p.flag_fortran += " -ffpe-trap=invalid"

    # GFortran can sets default values to a some unitialized memory, but it
    # would silent '-Wunitialized' so it is not enabled.
    #-finit-character=65 -finit-integer=42424242 -finit-real=snan


    #############
    # Libraries #
    #############
    p.library_list += " " + " ".join([
        # Undefined behavior sanitizer library.
        "ubsan",

        # Asks to 'malloc', defined in Glibc, to check the consistency
        # of dynamic memory by using the 'mcheck' functions.
        # Is equivalent do defines env variable MALLOC_CHECK_ to 3.
        # *** DISABLED since not thread safe
#        "mcheck",

        # Used to generate a backtrace in Talos.
        "backtrace",
    ])

    return p


def debug_mem(utils):
    """Debugging options for memory bugs at the cost of a slower speed.

    This profile is a layer on the 'debug' profile, see 'debug' profile for
    documentation.

    Note that by default 'undefined_behavior=log' and 'float_div_by_zero=log'
    are given to the 'debug' profile since the GCC sanitizer does not allow
    to generate core dump (so aborting the program is not as useful).

    It should be faster than the Valgrind Memcheck tool, with a complementary
    detection ability.
    """
    # When the address sanitizer is enabled, it is not possible to
    # generate core dumps, so aborting is not a very useful default.
    p = utils.profile
    p.preprocessor_defines += " TALOS_NO_COREDUMP"
    if "undefined_behavior" not in ARGUMENTS:
        ARGUMENTS["undefined_behavior"] = "log"
    if "float_div_by_zero" not in ARGUMENTS:
        ARGUMENTS["float_div_by_zero"] = "log"

    debug(utils)

    #########################
    # GCC address sanitizer #
    #########################

    # See https://code.google.com/p/address-sanitizer/
    # and notably https://code.google.com/p/address-sanitizer/wiki/Flags#Run-time_flags
    # => The average slowdown of the instrumented program is ~2x

    # This is not available in old version (YK).
    # p.flag_compiler += " -fsanitize=address"

    # It is implied by -fsanitize=address
    #p.flag_compiler += " -fsanitize=leak"

    p.flag_link += " -fsanitize=address"

    # Asan library needs the Ubsan one, so it must be put ahead in the list.
    p.library_list = "asan " +  p.library_list

    return p


def debug_thread(utils):
    """Debugging options for thread bugs at the cost of a slower speed.

    This profile is a layer on the 'debug' profile, see 'debug' profile for
    documentation.

    Its detection abilities are complementary with the Valgrind Helgrind tool.
    """
    p = debug(utils)

    ########################
    # GCC thread sanitizer #
    ########################

    # See https://code.google.com/p/thread-sanitizer/
    # and notably https://code.google.com/p/thread-sanitizer/wiki/Flags
    # with https://code.google.com/p/thread-sanitizer/wiki/CppManual
    # => memory usage may increase by 5-10x and execution time by 2-20x.
    # Can be be used if neither undefined nor thread sanitizers are enabled.
    # See https://code.google.com/p/address-sanitizer/wiki/LeakSanitizer
    #-fsanitize=thread
    p.flag_compiler +=  " -fsanitize=thread"

    p.flag_link += " -fsanitize=thread"
    p.library_list = "tsan " +  p.library_list

    return p

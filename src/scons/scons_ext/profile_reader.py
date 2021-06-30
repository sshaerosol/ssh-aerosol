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

from subprocess import Popen, PIPE
import os, re, types

from SCons.Errors import UserError

from . import FlagDict
import scons_ext.profile


def load_profile(utils, env, profile_name):
    """Imports the variables defined by the profile."""

    compiler, available_profile = list_profile(env)
    if profile_name not in available_profile:
        msg = ""
        if profile_name != '?':
            msg += "[ERROR] No '{0}' profile found for '{1}'\n"\
                      .format(profile_name, compiler)
        if not available_profile:
            msg += "[WARNING] No available profile.\n"
        else:
            msg += "List of available profiles for '{0}':\n".format(compiler)
            profile_list = available_profile.items()
            profile_list.sort()
            for p_name, p_func in profile_list:
                if p_func.__doc__:
                    p_doc = p_func.__doc__
                else:
                    p_doc = "(no documentation yet.)"
                msg += """
 '{0}'
 {1}
{2}""".format(p_name, '-' * (len(p_name) + 2), p_doc)
        profile_dir = os.path.dirname(scons_ext.profile.__file__)
        msg += """

------------------------------------------------------------------------
|  IMPORTANT: Profiles can be added and enriched, just have a look in: |
|  "{0}"
----------------------------------------------------------------------/
""".format(profile_dir)
        raise UserError(msg)

    available_profile[profile_name](utils)

def getstatusoutput(command):
    process = Popen(command, stdout=PIPE)
    out, _ = process.communicate()
    return (process.returncode, out)    
    
def list_profile(env):
    compiler = env.subst("$CC")

    # In the case MPI is used, gets the wrapped compiler.
    if compiler == "mpicc":
        s, o = getstatusoutput([compiler, "--showme", "-v"])        
        if s == 0:
            compiler = o.split()[0]
    elif compiler == "mpiicc":
        s, o = getstatusoutput([compiler, " -show"])
        if s == 0:
            compiler = o.split()[0]

    compiler_alias = [compiler]

    # Strips the version suffix.
    compiler_alternative = compiler.split('-')[0]
    if compiler_alternative != compiler:
        compiler_alias.append(compiler_alternative)

    profile_dict = {}
    for c in compiler_alias:
        for m_name in dir(scons_ext.profile):
            m = getattr(scons_ext.profile, m_name)
            if type(m) != types.ModuleType:
                continue
            if m_name != c:
                continue
            for p_name in dir(m):
                if p_name[0] == '_' or p_name in profile_dict:
                    continue
                p = getattr(m, p_name)
                try:
                    if p.__module__ != m.__name__:
                        continue
                except:
                    continue
                profile_dict[p_name] = p
    return (compiler, profile_dict)

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


class FlagDict:
    """Holds a profile attributes."""
    def __init__(self):
        self.__dict__['_dict'] = {}

    def __getattr__(self, name):
        """Returns an empty string when an attribute was not found."""
        try:
            return self._dict[name]
        except KeyError:
            default = ""
            self._dict[name] = default
            return default

    def __setattr__(self, name, value):
        self._dict[name] = value

    def dict(self):
        return self._dict

    def get(self, name, default=None):
        return self._dict.get(name, default)


# Example of use:
if __name__ == '__main__':
    example = FlagDict()
    example.a_flag += " -Wall"
    example.a_flag += " -Wpedantic"
    print(example.a_flag)
    print(example.dict())

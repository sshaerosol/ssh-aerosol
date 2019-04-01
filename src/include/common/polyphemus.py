#!/usr/bin/env python

# Copyright (C) 2006-2007, ENPC - INRIA - EDF R&D
#     Author(s): Vivien Mallet
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


# This module provides facilities to launch Polyphemus programs.


try:
    from network import Network
except:
    Network = str


##############
# POLYPHEMUS #
##############


class Polyphemus:
    """
    This class manages the execution of several Polyphemus programs.
    """


    def __init__(self, net = Network("comp")):
        """
        Initializes the network and the logs.

        @type net: Network instance
        @param net: The network over which the simuations should be launched.
        """
        self.program_list = []
        self.log = "-" * 78 + "\n\n"
        self.process = []
        self.net = net


    def GetLog(self):
        """
        Returns simulation logs.
        """
        return self.log


    def SetNetwork(self, net = Network("comp")):
        """
        Sets the network.

        @type net: Network instance
        @param net: The network over which the simuations should be launched.
        """
        self.net = net


    def AddProgram(self, program):
        """
        Adds a program.

        @type program: Program instance
        @param program: The program to be added.
        """
        if isinstance(program, str):
            self.program_list.append(Program(program))
        else:
            self.program_list.append(program)
        def compare(x, y):
            if x.group < y.group:
                return -1
            elif x.group > y.group:
                return 1
            else:
                return 0
        self.program_list.sort(compare)


    def Clear(self):
        """
        Clears all; primarily, the program list and the logs.
        """
        self.program_list = []
        self.log = "-" * 78 + "\n\n"


    def Run(self):
        """
        Executes the set of programs on localhost.
        """
        for program in self.program_list:
            program.Run()
            self.log += program.log
            if self.log[-1] != "\n":
                self.log += "\n"
            self.log += "\n" + "-" * 78 + "\n\n"
            if program.status != 0:
                print self.log
                raise Exception, "Program \"" + program.basename \
                      + "\" failed."


    def RunNetwork(self, Nhost = 10, delay = 30):
        """
        Executes the set of programs on the network.

        @type Nhost: integer
        @param Nhost: The number of hosts that may be used simultaneouly.
        @type delay: float or integer
        @param delay: The minimum period of time between the launch of two
        programs. Unit: seconds.
        """
        import time
        self.process = []
        self.host = []
        self.beg_time = []
        self.end_time = ["" for i in range(len(self.program_list))]
        # Index of the first program from current group.
        i_group = 0
        for i in range(len(self.program_list)):
            program = self.program_list[i]
            if i > i_group and program.group != self.program_list[i-1].group:
                # If any process from the previous group is still up.
                while min([x.poll() for x in self.process[i_group:]]) == -1:
                    time.sleep(delay)
                i_group = i
            host = self.net.GetAvailableHost(wtime = 2)
            p = self.net.LaunchBG(program.Command(), host = host)
            self.process.append(p)
            self.host.append(host)
            self.beg_time.append(time.asctime())

            while sum([1 for x in self.process[i_group:] if x.poll() == -1]) \
                      >= Nhost:
                # Maximum load reached.
                time.sleep(delay)

            for j in range(i):
                if self.end_time[j] == "" and self.process[i].poll() != 1:
                    self.end_time[j] = time.asctime()

        # Waits for the latest programs.
        while min([x.poll() for x in self.process[i_group:]]) == -1:
            time.sleep(delay)
            for j in range(len(self.program_list)):
                if self.end_time[j] == "" and self.process[i].poll() != 1:
                    self.end_time[j] = time.asctime()

        i_group = 0
        for i in range(len(self.program_list)):
            program = self.program_list[i]
            # New group ?
            if i > i_group and program.group != self.program_list[i-1].group:
                self.log += ("### GROUP " + str(program.group)
                             + " ###").center(78)
                self.log += "\n\n" + "-" * 78 + "\n\n"
                i_group = i

            self.log += program.Command()
            if self.log[-1] != "\n":
                self.log += "\n"
            self.log += "\nStatus: " + str(self.process[i].poll()) + "\n"
            self.log += "Hostname: " + str(self.host[i]) + "\n"
            self.log += "Started at " + str(self.beg_time[i]) + "\n"
            self.log += "Ended approximatively at " + str(self.end_time[i]) \
                        + "\n"
            self.log += "\n" + "-" * 78 + "\n\n"


    def Try(self):
        """
        Performs a dry run.
        """
        for program in self.program_list:
            program.Try()
            self.log += program.log
            if self.log[-1] != "\n":
                self.log += "\n"
            self.log += "\n" + "-" * 78 + "\n\n"
            if program.status != 0:
                raise Exception, "Program \"" + program.basename \
                      + "\" failed (status " + str(program.status) + ")."


###########
# PROGRAM #
###########


class Program:
    """
    This class manages a program associated with configuration files.
    """


    def __init__(self, name = None, config = None, format = " %a", group = 0):
        """
        Full initialization.

        @type name: string
        @param name: The program name.
        @type config: Configuration instance
        @param config: The program configuration.
        @type format: string
        @param format: The format of arguments, where "%a" is replaced with
        the configuration files.
        @type group: integer
        @param group: The group index.
        """
        if config is not None:
            self.config = config
        else:
            self.config = Configuration()
        self.name = name
        import os
        self.basename = os.path.basename(name)
        self.exec_path = "./"
        self.format = format
        self.priority = "0"
        self.output_directory = None
        self.group = group
        self.status = None
        self.log = None


    def Run(self):
        """
        Executes the program.
        """
        self.config.Proceed()
        import commands
        self.status, self.log = commands.getstatusoutput(self.Command())


    def Command(self):
        """
        Returns the command to launch the program. The program must be ready
        to be launched.
        """
        if not self.IsReady():
            raise Exception, "Program \"" + self.name + "\" is not ready."
        format = self.format[:]
        command = "nice time " + self.name \
                  + format.replace("%a", self.config.GetArgument())
        return command


    def SetConfiguration(self, config, mode = "random", path = None,
                         replacement = None, additional_file_list = []):
        """
        Sets the program configuration files.

        @type config: Configuration instance, or list of configuration files
        @param config: The configuration files associated with the program.
        @type mode: string
        @param mode: The copy mode. Mode "raw" just copies to the target path,
        while mode "random" appends a random string at the end of the file
        name. Mode "random_path" appends a random directory in the path. This
        entry is useless if "config" is a Configuration instance.
        @type path: string
        @param path: The path where configuration files should be copied. If
        set to None, then the temporary directory "/tmp" is used. This entry
        is useless if "config" is a Configuration instance.
        @type replacement: map with string keys
        @param replacement: The map of replaced strings and the replacement
        values. This entry is useless if "config" is a Configuration
        instance.
        @type additional_file_list: string or list of strings
        @param additional_file_list: An additional configuration file or a
        list of additional configuration files to be managed. Just like
        primary configuration files, they are subject to the replacements and
        copies, but are not considered as program arguments.
        """
        if isinstance(config, list) or isinstance(config, str):
            if mode is None:
                mode = "tmp"
            self.config = Configuration(config, mode, path,
                                        additional_file_list)
            if replacement is not None:
                self.config.SetReplacementMap(replacement)
        else:
            self.config = config


    def Try(self):
        """
        Performs a dry run.
        """
        self.config.Proceed()
        import os
        if os.path.isfile(self.name):
            self.status = 0
        else:
            self.status = 1
        format = self.format[:]
        command = self.name + format.replace("%a", self.config.GetArgument())
        self.log = "Running program \"" + self.basename + "\":\n" \
                   + "   " + command


    def IsReady(self):
        """
        Checks whether the program can be launched.

        @rtype: Boolean
        @return: True if the program can be executed, False otherwise.
        """
        return self.config.IsReady()


#################
# CONFIGURATION #
#################


class Configuration:
    """
    This class manages configuration files. It proceeds replacements in the
    files and makes copies of the files.
    """


    def __init__(self, file_list = [], mode = "random", path = None,
                 additional_file_list = []):
        """
        Initialization of configuration information.

        @type file_list: string or list of strings
        @param file_list: The configuration file or the list of configuration
        files to be managed.
        @type mode: string
        @param mode: The copy mode. Mode "raw" just copies to the target path,
        while mode "random" appends a random string at the end of the file
        name. Mode "random_path" appends a random directory in the path.
        @type path: string or None
        @param path: The path where configuration files should be copied. If
        set to None, then the temporary directory "/tmp" is used.
        @type additional_file_list: string or list of strings
        @param additional_file_list: An additional configuration file or a
        list of additional configuration files to be managed. Just like
        primary configuration files, they are subject to the replacements and
        copies, but are not considered as program arguments.
        """
        import os
        if isinstance(file_list, str):
            self.raw_file_list = [file_list]
        else:
            self.raw_file_list = file_list
        self.Narg = len(self.raw_file_list)
        if isinstance(additional_file_list, str):
            self.raw_file_list.append(additional_file_list)
        else:
            self.raw_file_list += additional_file_list
        for f in self.raw_file_list:
            if not os.path.isfile(f):
                raise Exception, "Unable to find \"" + f + "\"."
        self.file_list = []
        self.ready = False
        self.SetMode(mode)
        self.SetPath(path)
        self.config = {}


    def SetMode(self, mode = "random"):
        """
        @type mode: string
        @param mode: The copy mode. Mode "raw" just copies to the target path,
        while mode "random" appends a random string at the end of the file
        name. Mode "random_path" appends a random directory in the path.
        """
        if mode not in ["random", "random_path", "raw"]:
            raise Exception, "Mode \"" + str(mode) + "\" is not supported."
        self.mode = mode


    def SetPath(self, path):
        """
        @type path: string or None
        @param path: The path where configuration fiels should be copied. If
        set to None, then the temporary directory "/tmp" is used.
        """
        if path is not None:
            self.path = path
        else:
            self.path = "/tmp/"


    def IsReady(self):
        """
        Tests whether the configuration files are ready for use.

        @type: Boolean
        @return: True if the configuration files are ready for use, False
        otherwise.
        """
        return self.ready or self.raw_file_list == []


    def GetReplacementMap(self):
        """
        Returns the map of replaced strings and the replacement values.

        @rtype: map with string keys
        @return: The map of replaced strings and the replacement values.
        """
        return self.config


    def SetReplacementMap(self, config):
        """
        Sets the map of replaced strings and the replacement values.

        @type config: map with string keys
        @param config: The map of replaced strings and the replacement
        values.
        """
        self.config = config


    def SetConfiguration(self, config, mode = "random", path = None):
        """
        Iniitialization of configuration information, except file names.

        @type mode: string
        @param mode: The copy mode. Mode "raw" just copies to the target path,
        while mode "random" appends a random string at the end of the file
        name. Mode "random_path" appends a random directory in the path.
        @type path: string or None
        @param path: The path where configuration fiels should be copied. If
        set to None, then the temporary directory "/tmp" is used.
        """
        self.SetMode(mode)
        self.SetPath(path)
        self.SetReplacementMap(config)
        self.Proceed()


    def Proceed(self):
        """
        Proceeds replacement in configuration files and copy them.
        """
        import os, shutil, fileinput
        self.file_list = []
        if self.mode == "random_path" and self.raw_file_list is not []:
            import tempfile
            random_path = tempfile.mkdtemp(prefix = self.path)
        for f in self.raw_file_list:
            if self.mode == "raw":
                if os.path.dirname(f) == self.path:
                    raise Exception, "Error: attempt to overwrite" \
                          + " the raw configuration file \"" + f + "\"."
                name = os.path.join(self.path, os.path.basename(f))
                shutil.copy(f, name)
            elif self.mode == "random":
                import tempfile
                name = os.path.join(self.path, os.path.basename(f))
                fd, name = tempfile.mkstemp(prefix = name + "-")
                shutil.copy(f, name)
            elif self.mode == "random_path":
                name = os.path.join(random_path, os.path.basename(f))
                shutil.copy(f, name)
            self.file_list.append(name)

        if self.file_list != []:
            for line in fileinput.input(self.file_list, 1):
                new_line = line
                for i in self.config.keys():
                    new_line = new_line.replace(str(i), str(self.config[i]))
                if self.mode == "random_path":
                    new_line = new_line.replace("%random_path%", random_path)
                print new_line,
            fileinput.close()
        self.ready = True


    def GetRawFileList(self):
        """
        Returns the list of reference (or raw) configuration files.

        @rtype: list of strings
        @return: The list of reference (or raw) configuration files.
        """
        return self.raw_file_list


    def SetRawFileList(self, file_list):
        """
        Sets the list of reference (or raw) configuration files.

        @type file_list: list of strings
        @param file_list: The list of reference (or raw) configuration files.
        """
        self.raw_file_list = file_list
        self.ready = False


    def Clear(self):
        """
        Clears all, including configuration file names.
        """
        self.raw_file_list = []
        self.file_list = []
        self.ready = False
        self.mode = "random"
        self.path = "/tmp/"
        self.config = {}


    def GetArgument(self):
        """
        Returns the list of program arguments.

        @rtype: string
        @return: The list of program arguments aggregated in a string (and
        split by an empty space).
        """
        if self.IsReady():
            return " ".join(self.file_list[:self.Narg])
        else:
            raise Exception, "Not ready."


if __name__ == "__main__":
    simulation = Polyphemus(Network("loc"))

    program = Program("/bin/echo", format = " Hello World!")
    simulation.AddProgram(program)

    simulation.Run()

    print simulation.log

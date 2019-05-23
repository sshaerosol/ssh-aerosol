import distutils.sysconfig
env = Environment(SWIGFLAGS = ['-c++', '-python'],
                  CPPPATH = [distutils.sysconfig.get_python_inc()],
                  SHLIBPREFIX = "")
env.SharedLibrary('_talos.so', ['Talos.cpp', 'talos.i'])

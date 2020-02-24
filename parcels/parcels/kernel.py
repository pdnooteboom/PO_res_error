from parcels.codegenerator import KernelGenerator, LoopGenerator
from parcels.compiler import get_cache_dir
from parcels.tools.error import ErrorCode, recovery_map as recovery_base_map
from parcels.field import FieldOutOfBoundError, Field, SummedField, NestedField
from parcels.tools.loggers import logger
from parcels.kernels.advection import AdvectionRK4_3D
from os import path, remove
import numpy as np
import numpy.ctypeslib as npct
import time
from ctypes import c_int, c_float, c_double, c_void_p, byref
import _ctypes
from sys import platform, version_info
from ast import parse, FunctionDef, Module
import inspect
from copy import deepcopy
import re
from hashlib import md5
import math  # noqa
import random  # noqa


__all__ = ['Kernel']


re_indent = re.compile(r"^(\s+)")


def fix_indentation(string):
    """Fix indentation to allow in-lined kernel definitions"""
    lines = string.split('\n')
    indent = re_indent.match(lines[0])
    if indent:
        lines = [l.replace(indent.groups()[0], '', 1) for l in lines]
    return "\n".join(lines)


class Kernel(object):
    """Kernel object that encapsulates auto-generated code.

    :arg fieldset: FieldSet object providing the field information
    :arg ptype: PType object for the kernel particle

    Note: A Kernel is either created from a compiled <function ...> object
    or the necessary information (funcname, funccode, funcvars) is provided.
    The py_ast argument may be derived from the code string, but for
    concatenation, the merged AST plus the new header definition is required.
    """

    def __init__(self, fieldset, ptype, pyfunc=None, funcname=None,
                 funccode=None, py_ast=None, funcvars=None, c_include=""):
        self.fieldset = fieldset
        self.ptype = ptype
        self._lib = None

        # Derive meta information from pyfunc, if not given
        self.funcname = funcname or pyfunc.__name__
        if pyfunc is AdvectionRK4_3D:
            warning = False
            if isinstance(fieldset.W, Field) and fieldset.W.creation_log != 'from_nemo' and \
               fieldset.W._scaling_factor is not None and fieldset.W._scaling_factor > 0:
                warning = True
            if type(fieldset.W) in [SummedField, NestedField]:
                for f in fieldset.W:
                    if f.creation_log != 'from_nemo' and f._scaling_factor is not None and f._scaling_factor > 0:
                        warning = True
            if warning:
                logger.warning_once('Note that in AdvectionRK4_3D, vertical velocity is assumed positive towards increasing z.\n'
                                    '         If z increases downward and w is positive upward you can re-orient it downwards by setting fieldset.W.set_scaling_factor(-1.)')
        if funcvars is not None:
            self.funcvars = funcvars
        elif hasattr(pyfunc, '__code__'):
            self.funcvars = list(pyfunc.__code__.co_varnames)
        else:
            self.funcvars = None
        self.funccode = funccode or inspect.getsource(pyfunc.__code__)
        # Parse AST if it is not provided explicitly
        self.py_ast = py_ast or parse(fix_indentation(self.funccode)).body[0]
        if pyfunc is None:
            # Extract user context by inspecting the call stack
            stack = inspect.stack()
            try:
                user_ctx = stack[-1][0].f_globals
                user_ctx['math'] = globals()['math']
                user_ctx['random'] = globals()['random']
                user_ctx['ErrorCode'] = globals()['ErrorCode']
            except:
                logger.warning("Could not access user context when merging kernels")
                user_ctx = globals()
            finally:
                del stack  # Remove cyclic references
            # Compile and generate Python function from AST
            py_mod = Module(body=[self.py_ast])
            exec(compile(py_mod, "<ast>", "exec"), user_ctx)
            self.pyfunc = user_ctx[self.funcname]
        else:
            self.pyfunc = pyfunc

        if version_info[0] < 3:
            numkernelargs = len(inspect.getargspec(self.pyfunc).args)
        else:
            numkernelargs = len(inspect.getfullargspec(self.pyfunc).args)

        assert numkernelargs == 3, \
            'Since Parcels v2.0, kernels do only take 3 arguments: particle, fieldset, time !! AND !! Argument order in field interpolation is time, depth, lat, lon.'

        self.name = "%s%s" % (ptype.name, self.funcname)

        # Generate the kernel function and add the outer loop
        if self.ptype.uses_jit:
            kernelgen = KernelGenerator(fieldset, ptype)
            kernel_ccode = kernelgen.generate(deepcopy(self.py_ast),
                                              self.funcvars)
            self.field_args = kernelgen.field_args
            self.vector_field_args = kernelgen.vector_field_args
            fieldset = self.fieldset
            for f in self.vector_field_args.values():
                Wname = f.W.ccode_name if f.W else 'not_defined'
                for sF_name, sF_component in zip([f.U.ccode_name, f.V.ccode_name, Wname], ['U', 'V', 'W']):
                    if sF_name not in self.field_args:
                        if sF_name != 'not_defined':
                            self.field_args[sF_name] = getattr(f, sF_component)
            self.const_args = kernelgen.const_args
            loopgen = LoopGenerator(fieldset, ptype)
            if path.isfile(c_include):
                with open(c_include, 'r') as f:
                    c_include_str = f.read()
            else:
                c_include_str = c_include
            self.ccode = loopgen.generate(self.funcname, self.field_args, self.const_args,
                                          kernel_ccode, c_include_str)

            basename = path.join(get_cache_dir(), self._cache_key)
            self.src_file = "%s.c" % basename
            self.lib_file = "%s.%s" % (basename, 'dll' if platform == 'win32' else 'so')
            self.log_file = "%s.log" % basename

    def __del__(self):
        # Clean-up the in-memory dynamic linked libraries.
        # This is not really necessary, as these programs are not that large, but with the new random
        # naming scheme which is required on Windows OS'es to deal with updates to a Parcels' kernel.
        if self._lib is not None:
            _ctypes.FreeLibrary(self._lib._handle) if platform == 'win32' else _ctypes.dlclose(self._lib._handle)
            del self._lib
            self._lib = None
            if path.isfile(self.lib_file):
                [remove(s) for s in [self.src_file, self.lib_file, self.log_file]]

    @property
    def _cache_key(self):
        field_keys = "-".join(["%s:%s" % (name, field.units.__class__.__name__)
                               for name, field in self.field_args.items()])
        key = self.name + self.ptype._cache_key + field_keys + ('TIME:%f' % time.time())
        return md5(key.encode('utf-8')).hexdigest()

    def remove_lib(self):
        # Unload the currently loaded dynamic linked library to be secure
        if self._lib is not None:
            _ctypes.FreeLibrary(self._lib._handle) if platform == 'win32' else _ctypes.dlclose(self._lib._handle)
            del self._lib
            self._lib = None
        # If file already exists, pull new names. This is necessary on a Windows machine, because
        # Python's ctype does not deal in any sort of manner well with dynamic linked libraries on this OS.
        if path.isfile(self.lib_file):
            [remove(s) for s in [self.src_file, self.lib_file, self.log_file]]
            basename = path.join(get_cache_dir(), self._cache_key)
            self.src_file = "%s.c" % basename
            self.lib_file = "%s.%s" % (basename, 'dll' if platform == 'win32' else 'so')
            self.log_file = "%s.log" % basename

    def compile(self, compiler):
        """ Writes kernel code to file and compiles it."""
        with open(self.src_file, 'w') as f:
            f.write(self.ccode)
        compiler.compile(self.src_file, self.lib_file, self.log_file)
        logger.info("Compiled %s ==> %s" % (self.name, self.lib_file))

    def load_lib(self):
        self._lib = npct.load_library(self.lib_file, '.')
        self._function = self._lib.particle_loop

    def execute_jit(self, pset, endtime, dt):
        """Invokes JIT engine to perform the core update loop"""
        if len(pset.particles) > 0:
            assert pset.fieldset.gridset.size == len(pset.particles[0].xi), \
                'FieldSet has different amount of grids than Particle.xi. Have you added Fields after creating the ParticleSet?'
        for g in pset.fieldset.gridset.grids:
            g.cstruct = None  # This force to point newly the grids from Python to C
        # Make a copy of the transposed array to enforce
        # C-contiguous memory layout for JIT mode.
        for f in self.field_args.values():
            if not f.data.flags.c_contiguous:
                f.data = f.data.copy()
        for g in pset.fieldset.gridset.grids:
            if not g.depth.flags.c_contiguous:
                g.depth = g.depth.copy()
            if not g.lon.flags.c_contiguous:
                g.lon = g.lon.copy()
            if not g.lat.flags.c_contiguous:
                g.lat = g.lat.copy()
        fargs = [byref(f.ctypes_struct) for f in self.field_args.values()]
        fargs += [c_float(f) for f in self.const_args.values()]
        particle_data = pset._particle_data.ctypes.data_as(c_void_p)
        self._function(c_int(len(pset)), particle_data,
                       c_double(endtime), c_float(dt), *fargs)

    def execute_python(self, pset, endtime, dt):
        """Performs the core update loop via Python"""
        sign_dt = np.sign(dt)

        # back up variables in case of ErrorCode.Repeat
        p_var_back = {}

        for p in pset.particles:
            ptype = p.getPType()
            # Don't execute particles that aren't started yet
            sign_end_part = np.sign(endtime - p.time)
            if (sign_end_part != sign_dt) and (dt != 0):
                continue

            # Compute min/max dt for first timestep
            dt_pos = min(abs(p.dt), abs(endtime - p.time))
            while dt_pos > 1e-6 or dt == 0:
                for var in ptype.variables:
                    p_var_back[var.name] = getattr(p, var.name)
                try:
                    p.dt = sign_dt * dt_pos
                    res = self.pyfunc(p, pset.fieldset, p.time)
                except FieldOutOfBoundError as fse:
                    res = ErrorCode.ErrorOutOfBounds
                    p.exception = fse
                except Exception as e:
                    res = ErrorCode.Error
                    p.exception = e

                # Update particle state for explicit returns
                if res is not None:
                    p.state = res

                # Handle particle time and time loop
                if res is None or res == ErrorCode.Success:
                    # Update time and repeat
                    p.time += sign_dt * dt_pos
                    dt_pos = min(abs(p.dt), abs(endtime - p.time))
                    if dt == 0:
                        break
                    continue
                elif res == ErrorCode.Repeat:
                    # Try again without time update
                    for var in ptype.variables:
                        if var.name not in ['dt', 'state']:
                            setattr(p, var.name, p_var_back[var.name])
                    dt_pos = min(abs(p.dt), abs(endtime - p.time))
                    break
                else:
                    break  # Failure - stop time loop

    def execute(self, pset, endtime, dt, recovery=None, output_file=None):
        """Execute this Kernel over a ParticleSet for several timesteps"""

        def remove_deleted(pset):
            """Utility to remove all particles that signalled deletion"""
            indices = [i for i, p in enumerate(pset.particles)
                       if p.state in [ErrorCode.Delete]]
            if len(indices) > 0 and output_file is not None:
                output_file.write(pset[indices], endtime, deleted_only=True)
            pset.remove(indices)

        if recovery is None:
            recovery = {}
        recovery_map = recovery_base_map.copy()
        recovery_map.update(recovery)

        # Execute the kernel over the particle set
        if self.ptype.uses_jit:
            self.execute_jit(pset, endtime, dt)
        else:
            self.execute_python(pset, endtime, dt)

        # Remove all particles that signalled deletion
        remove_deleted(pset)

        # Identify particles that threw errors
        error_particles = [p for p in pset.particles
                           if p.state != ErrorCode.Success]
        while len(error_particles) > 0:
            # Apply recovery kernel
            for p in error_particles:
                if p.state == ErrorCode.Repeat:
                    p.state = ErrorCode.Success
                else:
                    recovery_kernel = recovery_map[p.state]
                    p.state = ErrorCode.Success
                    recovery_kernel(p, self.fieldset, p.time)

            # Remove all particles that signalled deletion
            remove_deleted(pset)

            # Execute core loop again to continue interrupted particles
            if self.ptype.uses_jit:
                self.execute_jit(pset, endtime, dt)
            else:
                self.execute_python(pset, endtime, dt)

            error_particles = [p for p in pset.particles
                               if p.state != ErrorCode.Success]

    def merge(self, kernel):
        funcname = self.funcname + kernel.funcname
        func_ast = FunctionDef(name=funcname, args=self.py_ast.args,
                               body=self.py_ast.body + kernel.py_ast.body,
                               decorator_list=[], lineno=1, col_offset=0)
        return Kernel(self.fieldset, self.ptype, pyfunc=None,
                      funcname=funcname, funccode=self.funccode + kernel.funccode,
                      py_ast=func_ast, funcvars=self.funcvars + kernel.funcvars)

    def __add__(self, kernel):
        if not isinstance(kernel, Kernel):
            kernel = Kernel(self.fieldset, self.ptype, pyfunc=kernel)
        return self.merge(kernel)

    def __radd__(self, kernel):
        if not isinstance(kernel, Kernel):
            kernel = Kernel(self.fieldset, self.ptype, pyfunc=kernel)
        return kernel.merge(self)

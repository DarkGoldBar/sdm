# -*- coding: utf-8 -*-
import sys
import time
import numpy as np
from enum import IntFlag
from typing import Tuple, Optional
from comlink.processor.base import SandboxOutput

from ..models import Lattice


class OptFlag(IntFlag):
    atom = 1
    volume = 2
    shape = 4


class OptimizeError(RuntimeError):
    pass


class FixSymmetry:
    def __init__(self, struct, space_group=None):
        self.struct = struct.getP1Structure()
        self.space_group = struct.space_group if space_group is None else space_group
        self.symm_map = None
        self.prep_symmetry()

    def prep_symmetry(self):
        fcoords = self.struct.frac_coords
        fcoords -= np.floor(fcoords)
        symm_map = []
        for op in self.space_group.symmetry_ops:
            this_op_map = [-1] * len(fcoords)
            for i_at in range(len(fcoords)):
                new_fc = op.rotation_matrix @ fcoords[i_at, :] + op.translation_vector
                dp = fcoords - new_fc
                dp -= np.floor(dp)
                i_at_map = np.argmin(np.sum(dp**2, axis=1))  # 简化 norm
                this_op_map[i_at] = i_at_map
            symm_map.append(this_op_map)
        self.symm_map = symm_map

    def adjust_cell(self, cell):
        # stress should definitely be symmetrized as a rank 2 tensor
        # UnitCellFilter uses deformation gradient as cell DOF with steps
        # dF = stress.F^-T quantity that should be symmetrized is therefore dF .
        # F^T assume prev F = I, so just symmetrize dF
        # F defined such that cell = cur_cell . F^T
        # assume prev F = I, so dF = F - I
        delta_deform_grad = np.dot(self.struct.lattice.inv_matrix, cell).T - np.eye(3)

        # symmetrization doesn't work properly with large steps, since
        # it depends on current cell, and cell is being changed by deformation
        # gradient
        max_delta_deform_grad = np.max(np.abs(delta_deform_grad))
        if max_delta_deform_grad > 0.25:
            raise RuntimeError('晶胞梯度过大 {} > 0.25'.format(max_delta_deform_grad))
        elif max_delta_deform_grad > 0.15:
            print('晶胞梯度较大 {} > 0.15'.format(max_delta_deform_grad))

        symmetrized_delta_deform_grad = symmetrize_rank2(self.struct.lattice,
                                                         delta_deform_grad,
                                                         self.space_group.symmetry_ops)
        cell[:] = self.struct.lattice.matrix @ (symmetrized_delta_deform_grad + np.eye(3)).T

    def adjust_coords(self, struct, new):
        step = new - struct.coords
        symmetrized_step = symmetrize_rank1(struct.lattice, step,
                                            self.space_group.symmetry_ops,
                                            self.symm_map)
        new[:] = struct.coords + symmetrized_step

    def adjust_forces(self, struct, forces):
        forces[:] = symmetrize_rank1(struct.lattice, forces,
                                     self.space_group.symmetry_ops,
                                     self.symm_map)


class BFGS:
    def __init__(self, struct, alpha=70.0, maxlstep=0.2, maxforce=0.001, maxnsteps=1000,
                 option: Optional[int] = None,
                 cell_pressure: Optional[float] = None,
                 cell_mask: Optional[Tuple[float]] = None,
                 fixsym: Optional[FixSymmetry] = None):
        self.struct = struct.getP1Structure() if hasattr(struct, 'space_group') else struct
        self.cell_mask = cell_mask
        self.cell_pressure = cell_pressure
        self.option = option
        self._drop_cache()

        # 优化设置
        if option is None:
            self.option = option = 7 if hasattr(struct, 'lattice') else 1
        if option > 1:
            if not hasattr(struct, 'lattice'):
                raise ValueError('输入对象没有晶胞 "lattice"!')
            self.cell_factor = len(struct.coords)
            self.cell_origin = struct.lattice.matrix.copy()
        elif option == 0:
            raise ValueError('错误的优化设置!')
        else:
            self.cell_factor = None
            self.cell_origin = None

        # 优化设置
        self.fixsym = fixsym
        self.alpha = alpha
        self.maxlstep = maxlstep
        self.maxnsteps = maxnsteps
        self.maxforce = maxforce

        # BFGS内部参数
        self.nsteps = 0
        self.H0 = np.eye(self.positions.size) * alpha
        self.H = None
        self.r0 = None
        self.f0 = None

    def _drop_cache(self):
        self._output = None
        self._forces = None
        self._positions = None
        self._cell_grad = None

    @property
    def output(self):
        if not self._output:
            self._output = self.calc(self.struct.toYoda())
        return self._output

    @property
    def cell_grad(self):
        if self._cell_grad is None:
            self._cell_grad = np.linalg.solve(self.cell_origin, self.struct.lattice.matrix).T
        return self._cell_grad

    @property
    def forces(self):
        if self._forces is None:
            if self.option > 1:
                stress = self.output.stress
                atoms_forces = self.output.forces @ self.cell_grad
                volume = self.struct.volume

                if self.cell_pressure:
                    stress = stress + np.eye(3) * self.cell_pressure

                virial = -volume * stress
                virial = np.linalg.solve(self.cell_grad, virial.T).T

                if not (self.option & OptFlag.shape):
                    vtr = virial.trace()
                    virial = np.eye(3) * (vtr / 3.0)

                if self.cell_mask is not None:
                    virial *= self.cell_mask

                if not (self.option & OptFlag.volume):
                    vtr = virial.trace()
                    np.fill_diagonal(virial, np.diag(virial) - vtr / 3.0)

                natoms = len(self.struct)
                forces = np.zeros((natoms + 3, 3))
                forces[:natoms] = atoms_forces
                forces[natoms:] = virial / self.cell_factor

                self.stress = virial / -volume
                self._forces = forces
            else:
                self._forces = self.output.forces
        return self._forces

    @property
    def positions(self):
        if self._positions is None:
            if self.option > 1:
                natoms = len(self.struct)
                pos = np.zeros((natoms + 3, 3))
                pos[:natoms] = np.linalg.solve(self.cell_grad, self.struct.coords.T).T
                pos[natoms:] = self.cell_factor * self.cell_grad
                self._positions = pos
            else:
                self._positions = self.struct.coords
        return self._positions

    @positions.setter
    def positions(self, val):
        if self.option > 1:
            natoms = len(self.struct)
            new_atom_positions = val[:natoms]
            new_deform_grad = val[natoms:] / self.cell_factor
            new_matrix = self.cell_origin @ new_deform_grad.T
            new_coords = new_atom_positions @ new_deform_grad.T
            if self.fixsym:
                self.fixsym.adjust_cell(self.struct, new_matrix)
                self.fixsym.adjust_coords(self.struct, new_coords)
            self.struct.lattice = Lattice(new_matrix)
            self.struct.coords = new_coords
        else:
            new_coords = val.copy()
            if self.fixsym:
                self.fixsym.adjust_coords(self.struct, new_coords)
            self.struct.coords = new_coords
        self._drop_cache()

    def calc(self) -> SandboxOutput:
        raise '运行时重载这个函数 self.calc(struct) -> SandboxOutput'

    def step(self, f=None):
        if f is None:
            f = self.forces
        r = self.positions
        f = f.reshape(-1)
        self.update(r.flat, f, self.r0, self.f0)

        omega, V = np.linalg.eigh(self.H)
        dr = np.dot(V, np.dot(f, V) / np.fabs(omega)).reshape((-1, 3))
        max_len_dr = np.sqrt(np.max((dr**2).sum(axis=1)))
        if max_len_dr >= self.maxlstep:
            scale = self.maxlstep / max_len_dr
            dr *= scale

        self.positions = r + dr
        self.r0 = r.flat.copy()
        self.f0 = f.copy()

    def update(self, r, f, r0, f0):
        if self.H is None:
            self.H = self.H0
            return
        dr = r - r0

        if np.abs(dr).max() < 1e-7:
            # Same configuration again (maybe a restart):
            return

        df = f - f0
        a = np.dot(dr, df)
        dg = np.dot(self.H, dr)
        b = np.dot(dr, dg)
        self.H -= np.outer(df, df) / a + np.outer(dg, dg) / b

    def log(self, forces=None, logfile=sys.stdout):
        if forces is None:
            forces = self.forces
        fmax = np.sqrt((forces ** 2).sum(axis=1).max())
        e = self.output.energy
        T = time.localtime()
        if logfile is not None:
            name = self.__class__.__name__
            if self.nsteps == 0:
                args = (" " * len(name), "Step", "Time", "Energy", "fmax")
                msg = "%s  %4s %8s %15s %12s\n" % args
                logfile.write(msg)

            args = (name, self.nsteps, T[3], T[4], T[5], e, fmax)
            msg = "%s:  %3d %02d:%02d:%02d %15.6f %12.4f\n" % args
            logfile.write(msg)

    def converged(self, forces=None):
        """Did the optimization converge?"""
        if forces is None:
            forces = self.forces
        return (forces ** 2).sum(axis=1).max() < self.maxforce ** 2

    def run(self):
        self._drop_cache()
        if self.nsteps == 0:
            self.log()
        while not self.converged() and self.nsteps < self.maxnsteps:
            self.step()
            self.nsteps += 1
            self.log()
        return self.converged()


def symmetrize_rank1(lattice, forces, symmetry_ops, symm_map):
    """
    Return symmetrized forces

    lattice vectors expected as row vectors (same as ASE get_cell() convention),
    inv_lattice is its matrix inverse (reciprocal().T)
    """
    scaled_symmetrized_forces_T = np.zeros(forces.T.shape)

    scaled_forces_T = np.dot(lattice.inv_matrix.T, forces.T)
    for (op, this_op_map) in zip(symmetry_ops, symm_map):
        transformed_forces_T = np.dot(op.rotation_matrix, scaled_forces_T)
        scaled_symmetrized_forces_T[:, this_op_map] += transformed_forces_T
    scaled_symmetrized_forces_T /= len(symmetry_ops)
    symmetrized_forces = (lattice.matrix.T @ scaled_symmetrized_forces_T).T

    return symmetrized_forces


def symmetrize_rank2(lattice, stress, symmetry_ops):
    """
    Return symmetrized stress

    lattice vectors expected as row vectors (same as ASE get_cell() convention),
    inv_lattice is its matrix inverse (reciprocal().T)
    """
    scaled_stress = np.dot(np.dot(lattice.matrix, stress), lattice.matrix.T)

    symmetrized_scaled_stress = np.zeros((3, 3))
    for op in symmetry_ops:
        symmetrized_scaled_stress += np.dot(np.dot(op.rotation_matrix.T, scaled_stress), op.rotation_matrix)
    symmetrized_scaled_stress /= len(symmetry_ops)

    sym = np.dot(np.dot(lattice.inv_matrix, symmetrized_scaled_stress), lattice.inv_matrix.T)
    return sym

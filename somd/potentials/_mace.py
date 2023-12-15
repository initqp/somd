#
# SOMD is an ab-initio molecular dynamics package designed for the SIESTA code.
# Copyright (C) 2023 github.com/initqp
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import os as _os
import numpy as _np
import typing as _tp
from somd import core as _mdcore
from somd.utils import constants as _c

__all__ = ['MACE']


class MACE(_mdcore.potential_base.POTENTIAL):
    """
    The MACE machine learning  potential [1].

    Parameters
    ----------
    atom_list : List[int]
        Indices of atoms included by this potential.
    file_name : str
        Name of the claculation method.
    atomic_types : List[int]
        Types of atoms included by this potential.
    device : str
        Name of the device to run MACE on.
    energy_unit : float
        The energy unit of the model outputs. In unit of (kJ/mol). E.g., if the
        model output energy in unit of (eV), this parameter should be set to
        96.485.
    length_unit : float
        The length unit of the model inputs and outputs. In unit of (nm). E.g.,
        if the model takes atomic positions in unit of (A) as input, this
        parameter should be set to 0.1.

    References
    ----------
    .. [1] Batatia, Ilyes, et al. "The design space of E(3)-equivariant
           atom-centered interatomic potentials." arXiv preprint arXiv (2022):
           2205.06643
    """

    def __init__(self,
                 atom_list: list,
                 file_name: str,
                 atomic_types: list,
                 device: str = 'cpu',
                 energy_unit: float = _c.AVOGACONST * _c.ELECTCONST * 0.001,
                 length_unit: float = 0.1) -> None:
        """
        Create a MACE instance.
        """
        super().__init__(atom_list)
        self.__energy_unit = energy_unit
        self.__length_unit = length_unit
        self.__atomic_types = _np.array(atomic_types, dtype=int).reshape(-1)
        self.__input_pbc = _np.array([True] * 3, dtype=bool)
        self.__input_forces = _np.zeros((len(atom_list), 3), dtype=float)
        self.__input_virial = _np.zeros((3, 3), dtype=float)
        self.__input_stress = _np.zeros(6, dtype=float)
        self.__input_charges = _np.zeros(len(atom_list), dtype=float)
        # treat MACE as a local dependency
        try:
            import torch
            import mace
            import mace.data
            import mace.tools
            self.__mace = mace
        except:
            raise ImportError('You need to have the mace package ' +
                              '(https://github.com/ACEsuit/mace) ' +
                              'installed to use the MACE potential!')
        model = torch.load(f=file_name, map_location=device)
        dtype = next(model.parameters()).dtype
        if (dtype == torch.float64):
            self.__model = model.to(device).double()
            self.__dtype = 'float64'
        elif (dtype == torch.float32):
            self.__model = model.to(device).float()
            self.__dtype = 'float32'
        else:
            message = 'Unknown model dtype: "{:s}"'.format(str(dtype))
            raise RuntimeError(message)
        mace.tools.torch_tools.set_default_dtype(self.__dtype)
        self.__r_max = float(model.r_max.cpu().numpy())
        self.__z_table = mace.tools.utils.AtomicNumberTable(
            [int(z) for z in model.atomic_numbers])
        for parameter in self.__model.parameters():
            parameter.requires_grad = False
        self.__device = mace.tools.torch_tools.init_device(device)

    def update(self, system: _mdcore.systems.MDSYSTEM) -> None:
        """
        Update this potential.

        Parameters
        ----------
        system : somd.systems.MDSYSTEM
            The simulated system.
        """
        configure = self.__mace.data.utils.Configuration(
            atomic_numbers=self.__atomic_types,
            positions=system.positions[self.atom_list] / self.__length_unit,
            energy=0.0,
            forces=self.__input_forces,
            stress=self.__input_stress,
            virials=self.__input_virial,
            dipole=None,
            charges=self.__input_charges,
            weight=1.0,
            energy_weight=0.0,
            forces_weight=0.0,
            stress_weight=0.0,
            virials_weight=0.0,
            config_type='Default',
            pbc=self.__input_pbc,
            cell=system.box / self.__length_unit)
        data_set = self.__mace.data.AtomicData.from_config(
            configure, z_table=self.__z_table, cutoff=self.__r_max)
        data_loader = self.__mace.tools.torch_geometric.dataloader.DataLoader(
            dataset=[data_set], batch_size=1, shuffle=False, drop_last=False)
        batch = next(iter(data_loader)).to(self.__device)
        result = self.__model(batch.to_dict(), compute_virials=True)
        energy = result['energy'].detach().cpu().numpy()
        self.energy_potential[0] = energy * self.__energy_unit
        forces = result['forces'].detach().cpu().numpy()
        self.forces[:] = forces * self.__energy_unit / self.__length_unit
        virial = result['virials']
        self.virial[:] = virial * self.__energy_unit

    @classmethod
    def generator(cls, *args, **kwargs) -> _tp.Callable:
        """
        Return a generator of this potential.
        """
        if 'file_name' in kwargs.keys():
            kwargs['file_name'] = _os.path.abspath(kwargs['file_name'])
        else:
            args = list(args)
            args[1] = _os.path.abspath(args[1])
        return lambda x=tuple(args), y=kwargs: cls(*x, **y)

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
from somd.utils import defaults as _d
from somd.utils import constants as _c
from somd.utils.warning import warn as _warn

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
    model_dtype : str
        Data type of the model.
    calculate_virial : bool
        If calculate the virial tensor.
    compile_mode : str
        Mode of compiling the model (see torch.compile for details).
    compile_full_graph : bool
        The `fullgraph` parameter of torch.compile.

    References
    ----------
    .. [1] Batatia, Ilyes, et al. "The design space of E(3)-equivariant
           atom-centered interatomic potentials." arXiv preprint arXiv (2022):
           2205.06643
    """

    def __init__(
        self,
        atom_list: _tp.List[int],
        file_name: str,
        atomic_types: _tp.List[int],
        device: str = 'cpu',
        energy_unit: float = _c.AVOGACONST * _c.ELECTCONST * 0.001,
        length_unit: float = 0.1,
        model_dtype: str = 'float64',
        calculate_virial: bool = True,
        compile_mode: str = None,
        compile_full_graph: bool = True,
    ) -> None:
        """
        Create a MACE instance.
        """
        super().__init__(atom_list)
        self.__energy_unit = energy_unit
        self.__length_unit = length_unit
        self.__calculate_virial = calculate_virial
        self.__comile_mode = compile_mode
        self.__full_graph = compile_full_graph
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
            self.__mace = mace  # fmt: skip
            self.__torch = torch
        except:
            raise ImportError(
                'You need to have the mace package '
                + '(https://github.com/ACEsuit/mace) '
                + 'installed to use the MACE potential!'
            )
        self.__file_name = file_name
        model = torch.load(f=file_name, map_location=device)
        # find head
        try:
            self.__heads = model.heads
        except AttributeError:
            self.__heads = ['Default']
        # try to compile
        if compile_mode is not None:
            if calculate_virial:
                message = 'Can not calculate virial using a compiled model!'
                raise RuntimeError(message)
            try:
                from mace.tools.compile import prepare
                from mace.tools.scripts_utils import extract_model
                model = torch.compile(
                    prepare(extract_model)(model=model, map_location=device),
                    mode=compile_mode,
                    fullgraph=compile_full_graph,
                )
            except Exception as e:
                self.__comile_mode = None
                message = 'Can not compile model, reason: {}'.format(e)
                _warn(message)
        defuault_dtype = next(model.parameters()).dtype
        if model_dtype is None:
            dtype = defuault_dtype
        elif model_dtype in ['float32', 'float']:
            dtype = torch.float32
        elif model_dtype in ['float64', 'double']:
            dtype = torch.float64
        else:
            message = 'Unknown model dtype: {}!'
            raise RuntimeError(message.format(model_dtype))
        if dtype != defuault_dtype:
            message = 'Model dtype has been convert from {} to {} as required.'
            _warn(message.format(defuault_dtype, dtype))
        if dtype == torch.float64:
            self.__model = model.to(device).double()
            self.__dtype = 'float64'
        elif dtype == torch.float32:
            self.__model = model.to(device).float()
            self.__dtype = 'float32'
        else:
            message = 'Unknown model dtype: "{:s}"'.format(str(dtype))
            raise RuntimeError(message)
        mace.tools.torch_tools.set_default_dtype(self.__dtype)
        self.__r_max = float(model.r_max.cpu().numpy())
        self.__z_table = mace.tools.utils.AtomicNumberTable(
            [int(z) for z in model.atomic_numbers]
        )
        for parameter in self.__model.parameters():
            parameter.requires_grad = False
        self.__device = mace.tools.torch_tools.init_device(device)
        self.__grad_outputs = [self.__torch.ones(1).to(self.__device)]

    def summary(self) -> str:
        """
        Show information about the potential.
        """
        result = '{}\n'.format(self.__class__.__name__)
        result += '┣━ n_atoms: {}\n'.format(self.n_atoms)
        result += '┣━ file_name: {}\n'.format(self.__file_name)
        result += '┣━ mace_version: {}\n'.format(self.__mace.__version__)
        result += '┣━ heads: {}\n'.format(self.__heads)
        result += '┣━ dtype: {}\n'.format(self.__dtype)
        result += '┣━ device: {}\n'.format(self.__device)
        result += '┣━ r_max: {}\n'.format(self.__r_max)
        result += '┣━ z_table: {}\n'.format(self.__z_table)
        result += '┣━ compile_mode: {}\n'.format(self.__comile_mode)
        result += '┣━ compile_full_graph: {}\n'.format(self.__full_graph)
        if _d.VERBOSE:
            result += '┣━ atom_list: {}\n'.format(self.atom_list)
        result += '┗━ END'
        return result

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
            cell=system.box / self.__length_unit,
        )
        data_set = self.__mace.data.AtomicData.from_config(
            configure,
            z_table=self.__z_table,
            cutoff=self.__r_max,
            heads=self.__heads
        )
        data_loader = self.__mace.tools.torch_geometric.dataloader.DataLoader(
            dataset=[data_set], batch_size=1, shuffle=False, drop_last=False
        )
        batch = next(iter(data_loader)).to(self.__device).to_dict()
        if self.__comile_mode is not None:
            batch['node_attrs'].requires_grad_(True)
            batch['positions'].requires_grad_(True)
        outputs = self.__model(
            batch,
            training=(self.__comile_mode is not None),
            compute_force=False,
            compute_virials=False,
            compute_stress=False,
            compute_displacement=self.__calculate_virial,
        )
        if self.__calculate_virial:
            forces, virial = self.__torch.autograd.grad(
                outputs=[outputs['energy']],
                inputs=[batch['positions'], outputs['displacement']],
                grad_outputs=self.__grad_outputs,
                retain_graph=False,
                create_graph=False,
                allow_unused=True,
            )
        else:
            forces = self.__torch.autograd.grad(
                outputs=[outputs['energy']],
                inputs=[batch['positions']],
                grad_outputs=self.__grad_outputs,
                retain_graph=False,
                create_graph=False,
                allow_unused=True,
            )
            forces = forces[0]
        if self.__calculate_virial and (virial is None):
            raise RuntimeError('Can not calculate virial!')
        if forces is None:
            raise RuntimeError('Can not calculate forces!')
        energy = outputs['energy'].detach().cpu().numpy()
        self.energy_potential[0] = energy * self.__energy_unit
        forces = forces.detach().cpu().numpy() * -1
        self.forces[:] = forces * self.__energy_unit / self.__length_unit
        if self.__calculate_virial:
            virial = virial.detach().cpu().numpy() * -1
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

import numpy as _np
import typing as _tp
from somd import core as _mdcore
from somd.utils import constants as _c

__all__ = ['DFTD4']


class DFTD4(_mdcore.potential_base.POTENTIAL):
    """
    The charge dependent London-dispersion correction potential [1].

    Parameters
    ----------
    atom_list : List[int]
        Indices of atoms included by this potential.
    atomic_types : List[int]
        Types of atoms included by this potential.
    method : str
        Name of the claculation method.
    total_charge : int
        Total charges of the atoms included by this potential.
    atm : bool
        If use the three-body correction.

    References
    ----------
    .. [1] Caldeweyher, Eike, et al. "A generally applicable atomic-charge
           dependent London dispersion correction." The Journal of chemical
           physics 150.15 (2019): 154122.
    """

    def __init__(
        self,
        atom_list: _tp.List[int],
        atomic_types: _tp.List[int],
        method: str,
        total_charge: int = 0,
        atm: bool = False,
    ) -> None:
        """
        Create a DFTD4 instance.
        """
        super().__init__(atom_list)
        # treat DFTD4 as a local dependency
        try:
            from dftd4.interface import DampingParam
            from dftd4.interface import DispersionModel
        except:
            raise ImportError(
                'You need to have the dftd4-python package '
                + 'installed to use the DFTD4 potential!'
            )
        self.__conversion = _c.HARTREE / _c.BOHRRADIUS * -1.0
        pbc = _np.ones(3, dtype=_np.int_)
        lattice = _np.zeros((3, 3), dtype=_np.double)
        atomic_types = _np.array(atomic_types, dtype=_np.int_)
        positions = _np.zeros((atomic_types.shape[0], 3), dtype=_np.double)
        self.__param = DampingParam(method=method, atm=atm)
        self.__model = DispersionModel(
            atomic_types, positions, total_charge, lattice, pbc
        )

    def update(self, system: _mdcore.systems.MDSYSTEM) -> None:
        """
        Update this potential.

        Parameters
        ----------
        system : somd.systems.MDSYSTEM
            The simulated system.
        """
        self.__model.update(
            system.positions[self.atom_list] / _c.BOHRRADIUS,
            system.box / _c.BOHRRADIUS,
        )
        result = self.__model.get_dispersion(self.__param, grad=True)
        self.energy_potential[0] = result.get("energy") * _c.HARTREE
        self.forces[:] = result.get("gradient") * self.__conversion
        self.virial[:] = result.get("virial") * _c.HARTREE * -1.0

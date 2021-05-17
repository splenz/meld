#
# Copyright 2015 by Justin MacCallum, Alberto Perez, Ken Dill
# All rights reserved
#

"""
A module to define the SystemState
"""

import numpy as np  # type: ignore


class SystemState:
    """
    Class to hold the state of a system.
    """

    def __init__(
        self,
        positions: np.ndarray,
        velocities: np.ndarray,
        alpha: float,
        energy: float,
        box_vector: np.ndarray,
    ) -> None:
        """
        Initialize a SystemState

        Params:
            positions: coordinates of structure, shape(n_atoms, 3)
            velocities: velocities for structure, shape(n_atoms, 3)
            alpha: alpha value, within ``[0, 1]``
            energy: total potential energy, including restraints, in kJ/mol
            box_vector: the box vectors, shape(3, 3) in nm
        """
        self.positions = positions
        self.velocities = velocities
        self.box_vector = box_vector
        self.n_atoms = positions.shape[0]
        self.alpha = alpha
        self.energy = energy

        self._validate()

    #
    # private methods
    #
    def _validate(self):
        # check positions
        if not len(self.positions.shape) == 2:
            raise RuntimeError("positions should be a 2D array")
        if not self.positions.shape[1] == 3:
            raise RuntimeError("positions should be (n_atoms, 3) array")

        # check velocities
        if not self.positions.shape == self.velocities.shape:
            raise RuntimeError("velocities must have the same shape as positions")

        # check box vectors
        if self.box_vector is not None:
            if not len(self.box_vector) == 3:
                raise RuntimeError("len(box_vectors) != 3")

        # check alpha
        if self.alpha < 0 or self.alpha > 1:
            raise RuntimeError("alpha must be in [0,1]")

import numpy as np

def fill_box_by_polymer(box, coords, atom_types=None, bonds=None, bond_types=None, angles=None, angle_types=None, 
                          dihedrals=None, dihedral_types=None, polymer_dist=1.0, dist_to_edge=0.1, nreplicate=None, full_box=True):
    '''
    Fill the box by replicating the polymer

    Parameters
    ----------
    box : list or np.ndarray
        Box size
    coords : np.ndarray
        Coordinates of the polymer
    atom_types : np.ndarray, optional
        Atom types of the polymer, by default None
    bonds : np.ndarray, optional
        Bonds of the polymer, by default None
    bond_types : np.ndarray, optional
        Bond types of the polymer, by default None
    angles : np.ndarray, optional
        Angles of the polymer, by default None
    angle_types : np.ndarray, optional
        Angle types of the polymer, by default None
    dihedrals : np.ndarray, optional
        Dihedrals of the polymer, by default None
    dihedral_types : np.ndarray, optional
        Dihedral types of the polymer, by default None
    polymer_dist : float, optional
        Distance between adjacent polymers, by default 1.0
    dist_to_edge : float, optional
        Distance between the polymer and the edge of the box, by default 0.1
    nreplicate : int, optional
        Number of replicate polymer, by default None
    full_box : bool, optional
        Whether to fill the box, by default True
    '''
    polymer_size = np.max(coords, axis=0) - np.min(coords, axis=0)
    polymer_atoms = coords.shape[0]
    box_size = np.array(box)
    if np.any(polymer_size>box_size):
        raise ValueError("Polymer size is larger than box size")
    coords = coords - box_size/2 + dist_to_edge - np.min(coords, axis=0) # move to the corner
    xyz_nrep = np.array([int((box_size[i]-2*dist_to_edge)//(polymer_size[i]+polymer_dist)) for i in range(3)]) # number of replicate polymer in each direction
    xyz_dist = (box_size-2*dist_to_edge)/xyz_nrep - polymer_size # distance between adjcent polymers
    if nreplicate is None and full_box:
        nreplicate = np.prod(xyz_nrep)
    else:
        if nreplicate > np.prod(xyz_nrep):
            raise ValueError("nreplicate is larger than the number of polymer can fill in the box")
        elif nreplicate < 1:
            raise ValueError("nreplicate should be larger than 1")
        elif not isinstance(nreplicate, int):
            raise ValueError("nreplicate should be an integer")

    shift_vecs = np.array([[i,j,k] for i in range(xyz_nrep[0]) 
                           for j in range(xyz_nrep[1]) for k in range(xyz_nrep[2])])[:nreplicate]*(polymer_size+xyz_dist)
    coords = np.tile(coords, (nreplicate, 1, 1)) + shift_vecs.reshape(-1, 1, 3)
    returns = [coords.reshape(-1, 3)]
    if atom_types is not None:
        atom_types = np.tile(atom_types, nreplicate)
        returns.append(atom_types)
    shift_atomid = (np.arange(nreplicate)*polymer_atoms).reshape(-1, 1, 1)
    if bonds is not None:
        bonds = np.tile(bonds, (nreplicate, 1, 1)) + shift_atomid
        returns.append(bonds.reshape(-1, 2))
    if bond_types is not None:
        bond_types = np.tile(bond_types, nreplicate)
        returns.append(bond_types)
    if angles is not None:
        angles = np.tile(angles, (nreplicate, 1, 1)) + shift_atomid
        returns.append(angles.reshape(-1, 3))
    if angle_types is not None:
        angle_types = np.tile(angle_types, nreplicate)
        returns.append(angle_types)
    if dihedrals is not None:
        dihedrals = np.tile(dihedrals, (nreplicate, 1, 1)) + shift_atomid
        returns.append(dihedrals.reshape(-1, 4))
    if dihedral_types is not None:
        dihedral_types = np.tile(dihedral_types, nreplicate)
        returns.append(dihedral_types)
    if len(returns)==1:
        return returns[0]
    else:
        return tuple(returns)

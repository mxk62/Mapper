from chemical import Chemical


def test_clear_ecs():
    """Returns True, if a molecule's atoms EC indices are not set."""

    # A chemical with no EC indices.
    chem = Chemical('CC1CCCC2CCCC(C)C12')
    chem.clear_ecs()
    assert chem.ec_order is None and not chem.ec_indices

    # A chemical with EC indices set.
    chem = Chemical('CCO')
    chem.ec_order = 1
    chem.ec_indices = (1, 2, 1)
    chem.clear_ecs()
    assert chem.ec_order is None and not chem.ec_indices


def test_find_ec_mcs():
    pass


def test_calc_init_ecs():
    """Returns True, if initial EC indices were calculated properly.

    Method calc_init_ecs() ONLY calculates initial EC indices, it does NOT
    alter the state of the class instance. Thus, attributes 'ec_order' and
    'ec_indices' should retain whatever values they had before calling it.
    """

    # A simple linear molecule, default Funatsu method.
    chem = Chemical('CCO')
    correct_indices = (61, 62, 81)
    assert chem.ec_order is None and not chem.ec_indices and \
           chem.calc_init_ecs() == correct_indices

    # A simple linear molecule, Shelley-Munk method.
    chem = Chemical('CCO')
    correct_indices = (21, 22, 41)
    assert chem.ec_order is None and not chem.ec_indices and \
           chem.calc_init_ecs(index_type='shelley') == correct_indices

    # A simple linear molecule, unknown method (defaults to Funatsu).
    chem = Chemical('CCO')
    correct_indices = (61, 62, 81)
    assert chem.ec_order is None and not chem.ec_indices and \
           chem.calc_init_ecs(index_type='undef') == correct_indices

    # Aliphatic fused rings, original Morgan method.
    chem = Chemical('CC1CCCC2CCCC(C)C12')
    correct_indices = (1, 3, 2, 2, 2, 3, 2, 2, 2, 3, 1, 3)
    assert chem.ec_order is None and not chem.ec_indices and \
           chem.calc_init_ecs(index_type='morgan') == correct_indices


def test_calc_next_ecs():
    """Returns True, if next order ECs were calculated properly.

    Method calc_next_ecs() ONLY calculates next order EC indices, it does NOT
    alter the state of the class instance. Thus, attributes 'ec_order' and
    'ec_indices' should retain whatever values they had before calling it.
    """

    # A simple linear molecule with' uninitialized indices.
    chem = Chemical('CCO')
    correct_indices = (184, 266, 224)
    assert chem.ec_order is None and not chem.ec_indices and \
           chem.calc_next_ecs() == correct_indices

    # A simple linear molecule with initialized with Morgan method.
    chem = Chemical('CCO')
    chem.ec_order = 1
    chem.ec_indices = (1, 2, 1)
    correct_indices = (4, 6, 4)
    assert chem.ec_order == 1 and chem.ec_indices == (1, 2, 1) and \
           chem.calc_next_ecs() == correct_indices


def test_remove_atoms():
    """Returns True, if atoms were removed properly."""

    # Removing terminating atom.
    chem = Chemical('CCO')
    chem.remove_atoms([0])
    expected_smiles = set(['CO'])
    generated_smiles = set([s for s in chem.smiles.split('.')])
    print generated_smiles
    assert generated_smiles.issubset(expected_smiles) == True

    # Removing atoms leaves disjoint fragments.
    chem = Chemical('CC(=O)CC(C)C(CC#N)C(=O)N')
    chem.remove_atoms([0, 1, 2, 3, 4, 5, 6])
    expected_smiles = set(['CC#N', 'NC=O'])
    generated_smiles = set([s for s in chem.smiles.split('.')])
    print generated_smiles
    assert generated_smiles.issubset(expected_smiles) == True

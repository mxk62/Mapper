from mapper import Chemical


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


def test_update_ecs():
    """Returns True, if EC indices were updated correctly."""

    # A chemical with having no EC indices.
    chem = Chemical('CCO')
    new_indices = (1, 2, 1)
    chem.update_ecs(new_indices)
    assert chem.ec_order == 1 and chem.ec_indices == new_indices

    # A chemical with initial EC indices set.
    chem = Chemical('CCO')
    chem.ec_order = 1
    chem.ec_indices = (1, 2, 1)
    new_indices = (4, 6, 4)
    chem.update_ecs(new_indices)
    assert chem.ec_order == 2 and chem.ec_indices == new_indices


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
    correct_indices = (12, 22, 14)
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

    # A simple linear molecule with' uninitialized indices (will default to
    # Funatsu method).
    chem = Chemical('CCO')
    correct_indices = (306, 390, 386)
    assert chem.ec_order is None and not chem.ec_indices and \
        chem.calc_next_ecs() == correct_indices

    # A simple linear molecule with initialized with Morgan method.
    chem = Chemical('CCO', ec_type='morgan')
    chem.ec_order = 1
    chem.ec_indices = (1, 2, 1)
    correct_indices = (4, 6, 4)
    assert chem.ec_order == 1 and chem.ec_indices == (1, 2, 1) and \
        chem.calc_next_ecs() == correct_indices

    # A simple linear molecule with initialized with Shelley-Munk method.
    chem = Chemical('CCO', ec_type='shelley')
    chem.ec_order = 1
    chem.ec_indices = (12, 22, 14)
    correct_indices = (46, 70, 50)
    assert chem.ec_order == 1 and chem.ec_indices == (12, 22, 14) and \
        chem.calc_next_ecs() == correct_indices


def test_remove_atoms():
    """Returns True, if atoms were removed properly."""

    # Removing terminating atom.
    chem = Chemical('CCO')
    chem.remove_atoms([0])
    expected_smiles = {'CO'}
    generated_smiles = set([s for s in chem.smiles.split('.')])
    assert generated_smiles.issubset(expected_smiles) is True

    # Removing atoms leaves disjoint fragments.
    chem = Chemical('CC(=O)CC(C)C(CC#N)C(=O)N')
    chem.remove_atoms([0, 1, 2, 3, 4, 5, 6])
    expected_smiles = {'CC#N', 'NC=O'}
    generated_smiles = set([s for s in chem.smiles.split('.')])
    assert generated_smiles.issubset(expected_smiles) is True

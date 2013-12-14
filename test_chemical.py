from chemical import Chemical


def test_clear_ecs():
    """Returns True, if a molecule's atoms has no 'EC' property."""

    # A chemical with no EC indices.
    chem = Chemical('CC1CCCC2CCCC(C)C12')
    chem.clear_ecs()
    assert any(a.HasProp('EC') for a in chem.mol.GetAtoms()) == False and \
        chem.ec_order is None

    # A chemical with EC indices set.
    chem = Chemical('CCO')
    for i, a in enumerate(chem.mol.GetAtoms()):
        a.SetProp('EC', '1' if i % 2 == 0 else '2')
    chem.clear_ecs()
    assert any(a.HasProp('EC') for a in chem.mol.GetAtoms()) == False and \
        chem.ec_order is None


def test_find_ecs_mcs():
    pass


def test_find_ecn():
    """Returns True, if all atoms EC indices were calculated property."""

    # A chemical with no EC indices.
    chem = Chemical('CCO')
    assert [chem.find_ecn(a) for a in chem.mol.GetAtoms()] == [4, 6, 4] \
        and chem.ec_order == 0

    # A chemical with EC indices set.
    chem = Chemical('CCO')
    chem.ec_order = 1
    for i, a in enumerate(chem.mol.GetAtoms()):
        a.SetProp('EC', '4' if i % 2 == 0 else '6')
    print chem.ec_order, [chem.find_ecn(a) for a in chem.mol.GetAtoms()]
    assert [chem.find_ecn(a) for a in chem.mol.GetAtoms()] == [14, 20, 14]


def test_init_ecs():
    """Returns True, if ECs were initialized properly."""

    # A simple linear molecule.
    chem = Chemical('CCO')
    assert chem.init_ecs() == (1, 2, 1) and chem.ec_order == 0

    # Aliphatic fused rings.
    chem = Chemical('CC1CCCC2CCCC(C)C12')
    assert chem.init_ecs() == (1, 3, 2, 2, 2, 3, 2, 2, 2, 3, 1, 3) and \
        chem.ec_order == 0


def test_increment_ecs():
    """Returns True, if ECs were incremented properly."""
    chem = Chemical('CCO')
    chem.ec_order = 0
    for i, a in enumerate(chem.mol.GetAtoms()):
        a.SetProp('EC', '1' if i % 2 == 0 else '2')
    assert chem.increment_ecs() == (4, 6, 4) and chem.ec_order == 1


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


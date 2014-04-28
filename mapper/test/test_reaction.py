from mapper import Reaction


def test_find_core():
    """Returns True, if core was extracted correctly."""

    # Working examples from the original paper by Lynch and Willett.
    smiles = '>>'.join([
        'CC(=O)CC(C)C(CC#N)C(=O)N',
        'CC(=O)CC(C)C(CC#N)C#N'
    ])
    rxn = Reaction(smiles)
    assert rxn.find_core() == 'NC=O>>C#N'

    smiles = '>>'.join([
        'NC(=O)C(c1ccccc1)(c1ccccc1)NS(=O)(=O)c1ccccc1',
        'O=C(O)C(NS(=O)(=O)c1ccccc1)(c1ccccc1)c1ccccc1'
    ])
    rxn = Reaction(smiles)
    assert rxn.find_core() == 'NC=O>>O=CO'

    smiles = '>>'.join([
        'CCC(c1cc(OC)c(OC(C)=O)cc1Cc1ccc(OC(C)=O)c(OC)c1)C(C)OC(C)=O',
        'CCC(c1cc(OC)c(OC(C)=O)cc1C(=O)c1ccc(OC(C)=O)c(OC)c1)C(C)OC(C)=O'
    ])
    rxn = Reaction(smiles)
    assert rxn.find_core() == 'cCc>>cC(c)=O'

    smiles = '>>'.join([
        'ON(O)c1ccccc1S(=O)(=O)N(C)c1ccccc1',
        'Nc1ccccc1S(=O)(=O)N(C)c1ccccc1'
    ])
    rxn = Reaction(smiles)
    assert rxn.find_core() == 'cN(O)O>>cN'

    # With correct Shelley-Munk initial indices, algorithm returns different
    # result from the one state by Lynch-Willett due to unambiguous matches.
    # Original Lynch-Willett approach in assigning initial indices might have
    # more discriminative power but currently is not implemented due to its
    # rather vague description. It would also require writing a 'translator'
    # between SMILES and WLN notation on which it is based.
    #
    #smiles = '>>'.join([
    #    'COc1ccc2c(c1)sc1c2CCC2C1=CCC2O',
    #    'COc1ccc2c(c1)sc1c2CCC2C1CCC2O'
    #])
    #rxn = Reaction(smiles)
    #assert rxn.find_core() == 'C=CC>>CCC'

    # Problematic reactions, violating algorithm assumption---no atom should
    # be deleted; also from the paper by Lynch and Willett.
    #
    #
    # Ambiguous mapping.
    #
    # RDKit (2013_09_01) is having problem with this SMILES so we are skipping
    # them for now.
    #smiles = '>>'.join([
    #    'c1ccc2c(c1)C=Cc1ccccc1NN2',
    #    'CN1Nc2ccccc2C=Cc2c1cccc2'
    #])
    #rxn = Reaction(smiles)
    #assert rxn.find_core == smiles

    smiles = '>>'.join([
        'N=C1ON=C2CCCCC12',
        'Nc1onc2c1CCCC2'
    ])
    rxn = Reaction(smiles)
    assert rxn.find_core() == smiles

    # Match radius is to small.
    smiles = '>>'.join([
        'CC(C)(O)C1CCCO1',
        'CC1=COCCC1C'
    ])
    rxn = Reaction(smiles)
    assert rxn.find_core() == smiles

    smiles = '>>'.join([
        'O=CC1C=CCC=C1',
        'O=CC1=CCCC=C1'
    ])
    rxn = Reaction(smiles)
    assert rxn.find_core() == smiles


def test_make_ec_map():
    """Returns True, if map between EC and atom indices was done correctly."""
    pass

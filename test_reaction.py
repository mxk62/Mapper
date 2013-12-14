from reaction import Reaction


def test_find_core():
    """Returns True, if core was extracted correctly."""

    # Test if it work properly at all.
    smiles = 'CC(=O)CC(C)C(CC#N)C(=O)N>>CC(=O)CC(C)C(CC#N)C#N'
    rxn = Reaction(smiles)
    assert rxn.find_core() == 'N#CC.C(=O)N>>N#CC.N#C'


def test_make_ec_map():
    """Returns True, if map between EC and atom indices was done correctly."""
    pass

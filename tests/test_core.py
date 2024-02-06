from mapper import Core


def test_does_contain():
    """Returns True if cores' inclusion in properly recognized."""

    # Single retron cores; the first core contain the other one.
    core = Core('[C:1][O:2]>>[C:1]=[O:2]')
    other = Core('[C:2]>>[C:2]Br')
    assert core.does_contain(other) is True

    # Single retron cores; the first core does not contain the other one.
    core = Core('[C:1][O:2]>>[C:1]=[O:2]')
    other = Core('[S:2]>>O=[S:2]=O')
    assert core.does_contain(other) is False

    # Multi-retron cores; the same numbers of retrons.
    core = Core('[C:4]O.[N:3][C:1]=[O:2]>>[C:4][O:2][C:1]=[N:3]')
    other = Core('[N:3].[C:2]O>>[C:2][N:3]')
    assert core.does_contain(other) is True

    # Multi-retron cores; different numbers of retrons.
    core = Core('[C:4]O.[N:3][C:1]=[O:2]>>[C:4][O:2][C:1]=[N:3]')
    other = Core('[C:1]=O>>[C:1]=C')
    assert core.does_contain(other) is False

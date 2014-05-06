from mapper import Pattern


def test_does_contain():
    """Returns True if cores' inclusion in properly recognized."""

    # Single fragment pattern; the first pattern contains the other one.
    patt = Pattern('[C:1][O:2]')
    other = Pattern('[C:2]')
    assert patt.does_contain(other) is True

    # Single fragment pattern; the first pattern does not contain the other.
    patt = Pattern('[C:1][O:2]')
    other = Pattern('[S:2]')
    assert patt.does_contain(other) is False

    # Multi-fragment patterns; the same numbers of fragments.
    patt = Pattern('[C:4]O.[N:3][C:1]=[O:2]')
    other = Pattern('[N:3].[C:2]O')
    assert patt.does_contain(other) is True

    # Multi-fragment patterns; different numbers of fragments.
    patt = Pattern('[C:4]O.[N:3][C:1]=[O:2]')
    other = Pattern('[C:1]=O')
    assert patt.does_contain(other) is False


def test_find_distance():
    """Returns True if similarity between patterns was evaluated correctly."""

    # Identical single-fragment patterns, distance 0.
    patt = Pattern('[C:1][O:2]')
    other = Pattern('[C:1][O:2]')
    assert abs(patt.find_similarity(other) - 0.0) < 1.0e-06

    # Partially overlapping, single-fragment patterns, distance 0.5.
    patt = Pattern('[C:1][O:2]')
    other = Pattern('[C:1]')
    assert abs(patt.find_similarity(other) - 0.5) < 1.0e-06

    # Non-overlapping, single-fragment patterns, distance 1.
    patt = Pattern('[C:1][O:2]')
    other = Pattern('[S:1]')
    assert abs(patt.find_similarity(other) - 1.0) < 1.0e-06

    # Identical multi-fragment patterns, distance 0.
    patt = Pattern('[C:4]O.[N:3][C:1]=[O:2]')
    other = Pattern('[C:4]O.[N:3][C:1]=[O:2]')
    assert abs(patt.find_similarity(other) - 0.0) < 1.0e-06

    # Partially overlapping, multi-fragment patterns, distance 0.5.
    patt = Pattern('[C:4]O.[N:3][C:1]=[O:2]')
    other = Pattern('[N:3].[C:2]O')
    assert abs(patt.find_similarity(other) - 0.333333) < 1.0e-06

    # Non-overlapping, multi-fragment patterns, distance 1.
    patt = Pattern('[C:4]O.[N:3][C:1]=[O:2]')
    other = Pattern('[S:1].[S:2]')
    assert abs(patt.find_similarity(other) - 1.0) < 1.0e-06

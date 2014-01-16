# mapper

Automatic reaction mapping and center detection


## Description

`mapper` is an implementation of approximate structure-matching algorithm
proposed by Lynch and Willett [1]. It applies the Morgan algorithm [2] with
different stopping rule *simultaneously* to both reacting molecules to detect
inter-molecular equivalences. 

The authors of Ref. [1] are a bit vague as far as assigning initial EC indices
is concerned. Thus, besides original Morgan's approach, two alternative methods
[3,4] were implemented.

## Manifest

Package source directory should contain:

 *  `README.md`: this file.
 *  `chemical.py`: definition of Chemical class;
 *  `reaction.py`: definition of Reaction class;
 *  `test_chemical.py`: unit tests for Chemical class;
 *  `test_reaction.py`: unit tests for Reaction class;
 * 	`demo.py`: a (very) simple CLI for demonstration purposes.


## References

1.	MF Lynch and P Willett, *The Automatic Detection of Chemical Reaction
	Sites*, J Chem Inf Comput Sci **18**, 154-159 (1978)
2.	HL Morgan, , *The generation of a Unique Machine Description for Chemical
	Structures---A Techinique Developed at Chemical Abstract Service*, J Chem
	Doc **5**, 107-113 (1965)
3.	CA Shelley and ME Munk, *Computer Perception of Topological Symmetry*,
	J Chem Inf Mod **17**, 110-113 (1977)
4.	K Funatsu, T Endo, N Kotera, and SI Sasaki, *Automatic Recognition of
	Reaction Site in Organic Chemnical Reactions*, Tetrahedron Computer
	Methodology **1**, 53-69 (1988)

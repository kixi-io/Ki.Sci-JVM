# Ki.Sci-JVM

The `sci` type bag: the models of matter that the physics, chemistry and biology
bags build on, for KD, KS, KQ and RevDB. A type bag is a library that adds value
types, functions and constants to every Ki layer at once (the Ki type bags design,
in Ki.Docs). Ki.Sci depends on Ki.Core alone.

```kd
compound name="ethanol" structure=.smiles("CCO") formula=.formula(C2H6O) inchi=.inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
```

```ks
use ki.sci.*
let m = .smiles("OCC")
say m == .smiles("CCO")       // true: compared by canonical SMILES
say m.formula                 // C2H6O
say m.mass                    // 46.069Da
say mass("NaCl")              // 58.43976928Da
say .element(Fe).group        // 8
```

```kq
from c in compounds where c.structure.formula == .formula(C2H6O) select c.name
```

## Types

| Type | Literal | Orders by |
|---|---|---|
| `Element` | `.element(Na)` (symbol, name or number) | atomic number |
| `Isotope` | `.isotope(13C)` (also `C-13`, `D`, `T`) | atomic number, then mass number |
| `Formula` | `.formula(C6H12O6)` (groups, hydrates, charges) | Hill text |
| `Molecule` | `.smiles("CCO")` | canonical SMILES |
| `InChI` | `.inchi("InChI=1S/…")` | text |

- **Elements**: the 118, with CIAAW 2021 standard atomic weights, their
  uncertainties and intervals, IUPAC's abridged weights, group, period and block.
- **Isotopes**: every nuclide in AME2020 (2,939), with natural abundances from
  CIAAW 2021. The data are in `src/main/resources/io/kixi/sci/isotopes.txt`.
- **Formulas**: Hill order, molar mass (standard weights, electrons not counted)
  and monoisotopic mass (main isotopes, electrons counted), composition.
- **Molecules**: SMILES (OpenSMILES; not wildcards, non-tetrahedral chirality or
  cumulene cis/trans) read to a graph with tetrahedral and double-bond
  stereochemistry, written back as canonical SMILES, which is what equality
  compares. Explicit hydrogens fold into counts (`[H]/C(F)=C/F` is `F/C=C\F`),
  stereo marks on atoms that are not stereocenters are dropped (`C[C@H](C)C` is
  `CC(C)C`), and an unwritten bond between aromatic rings is single
  (`c1ccccc1c1ccccc1`), as RDKit, Open Babel and CDK read them. There is no
  aromaticity model: aromatic and Kekulé forms of one molecule are different
  values until Ki.Chem.Ext (CDK) adds one. Up to 5,000 atoms.
- **InChI**: checked and kept as text. Keys need the InChI library, which
  Ki.Chem.Ext brings.

Functions: `element(x)`, `isotope(x)` or `isotope(element, massNumber)`,
`formula(x)`, `molecule(x)`, `mass(x)` (text is a formula: `mass("CCO")` is
C2O's mass; write `mass(molecule("CCO"))` for ethanol). Constants (CODATA
2022, as quantities): `speedOfLight`, `electronMass`, `protonMass`,
`neutronMass`. The constants
whose units Ki.Core does not have yet (Avogadro, Boltzmann, Planck, the gas
constant, Faraday, the elementary charge) wait for those units.

## Build

A composite build with the sibling checkouts `../Ki.Core-JVM`, `../Ki.KD-JVM`,
`../Ki.KS-JVM` and `../RevDB-JVM` (the last three for the tests only). Kotlin
2.3.0, Kotest 6.0.7, JDK 21.

```
./gradlew test
```

The canonical SMILES tests write every molecule of a corpus, and of two random
fuzz sets, in many random atom orders and check that each reads back to the same
canonical form. The corpus round trip puts every sample value through KD, KS
and RevDB.

MIT license. The element and isotope data are from IUPAC/CIAAW and the AME2020
evaluation, by way of the public-domain `periodictable` package.

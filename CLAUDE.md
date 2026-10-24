# Ki.Sci-JVM

The `sci` type bag (Ki type bags design, sections 2 and 5.2; the design is in Ki.Docs and the RevDB project docs). Kotlin 2.3.0 with `explicitApi()`, Kotest 6.0.7, JDK 21.

## Composite build
Includes `../Ki.Core-JVM` (the only main dependency) and `../Ki.KD-JVM`, `../Ki.KS-JVM`, `../RevDB-JVM` (tests only). Ask before adding a dependency.

## Rules
- The bag is `SciBag` (a class with a public no-argument constructor, listed in `META-INF/services/io.kixi.types.TypeBag`); nothing registers it. Tests reach it with `TypeBags.byName("sci")`.
- Stored encodings and sort keys never change for a released version: bump `SciBag.VERSION` and keep decoding older versions. A Molecule is stored and indexed as its canonical SMILES, so **any change to `SmilesReader`, `SmilesCanon` or `SmilesWriter` that changes a canonical string changes stored sort keys**: it needs a version bump and an index rebuild (an open question in the design, 9.4).
- After touching the SMILES code, run `MoleculeTest`: the corpus and the two fuzz sets write each molecule in many random atom orders and check one canonical form. Add a failing SMILES to the corpus.
- Element and isotope data come from CIAAW 2021 and AME2020 (`isotopes.txt` is generated; see its header). Change them only from those sources, never by hand.
- Units are not part of a bag (design 9.2, item 6): a constant whose unit Ki.Core lacks waits for Ki.Core.
- Function and constant names are global across bags (and KQ's built-ins win in KQ): check the other bags before adding one.
- `BagCorpus.kt` (the corpus round trip) is shared with Ki.Math-JVM's and Ki.Bio-JVM's tests; keep the copies alike.

## Working agreements
- Commit straight to `master`; never create a branch unless Dan asks. Code and docs in separate commits, author `dleuck <dan@kixi.io>`. Don't publish or change the version without asking.

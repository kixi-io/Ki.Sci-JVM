package io.kixi.sci

import io.kotest.assertions.throwables.shouldThrow
import io.kotest.assertions.withClue
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.shouldBe
import io.kotest.matchers.shouldNotBe
import io.kotest.matchers.string.shouldContain
import java.math.BigDecimal
import kotlin.random.Random

/** SMILES reading, canonical form (invariance under atom order, stereochemistry), and the derived values. */
class MoleculeTest : FunSpec({

    fun m(s: String) = Molecule.parse(s)
    fun near(q: io.kixi.uom.Quantity<*>, lo: String, hi: String) {
        val v = q.value as BigDecimal
        withClue("$q in [$lo, $hi]") { (v >= BigDecimal(lo) && v <= BigDecimal(hi)) shouldBe true }
    }
    fun canon(s: String) = m(s).smiles

    /** The molecule [s] written with random ranks: a SMILES of the same molecule in another atom order. */
    fun shuffled(s: String, random: Random): String {
        val g = SmilesReader.read(s)
        val ranks = (0 until g.n).shuffled(random).toIntArray()
        return SmilesWriter(g, ranks).write().smiles
    }

    val corpus = listOf(
        "C", "CC", "CCO", "OCC", "C(O)C", "CC(=O)O", "CC(C)(C)C", "C1CCCCC1", "c1ccccc1", "C1=CC=CC=C1",
        "c1ccc2ccccc2c1", "c1ccc2cc3ccccc3cc2c1", "C1CC2CCC1CC2", "C12C3C4C1C5C2C3C45", "O=C1C=CC(=O)C=C1",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "CC(=O)OC1=CC=CC=C1C(=O)O", "[Na+].[Cl-]", "O.O.O", "[NH4+].[NH4+].[O-]S(=O)(=O)[O-]",
        "[13CH4]", "[2H]O[2H]", "C[2H]", "[CH3:1][OH:2]", "N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O", "C[C@@H](C(=O)O)N",
        "F/C=C/F", "F/C=C\\F", "C/C=C/C=C/C", "C/C=C\\C=C/C", "CC/C=C(\\C)/C=C/C", "C[S@](=O)CC", "C[S@@](=O)CC",
        "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "C[C@]12CC[C@H]3[C@@H](CC=C4C[C@@H](O)CC[C@]34C)[C@@H]1CC[C@@H]2O",
        "F[C@H](Cl)Br", "F[C@@H](Cl)Br", "[C@@H](F)(Cl)Br", "Cl[C@H](F)Br", "F[C@H](F)Cl", "c1cc[nH]c1", "[O-][N+](=O)c1ccccc1",
        "C#N", "[Cu+2]", "[Fe+3].[Fe+3].[O-2].[O-2].[O-2]", "C1CC1C1CC1", "C1C/C=C/CCCC1", "c1ccccc1-c1ccccc1",
        "CCCCCCCCCCCCCCCCCCCC", "C(C(C(C(C(C)C)C)C)C)C", "O=C=O", "[H][H]", "[H+]", "[OH-]", "B(F)(F)F", "P(=O)(O)(O)O",
        "S(=O)(=O)(O)O", "ClC(Cl)(Cl)Cl", "C1=CC=C2C(=C1)C=CC=C2", "N#Cc1ccccc1", "C=C=C", "CC(C)[C@@H](C(=O)O)N",
        "C1[C@H]2CC[C@@H]1CC2", "[C@H]1(F)CC[C@@H](F)CC1", "F/C(Cl)=C(/F)Cl", "C(/F)=C/F", "C=1CCCC=1",
        // Found by the review: stereo alkenes on one carbonyl carbon, biphenyl without '-', a marked [H].
        "O=C(/C=C/C)/C=C/C", "C\\C=C/C(=O)/C=C\\C", "C/C=C\\C(=O)/C=C/C", "C\\C=C/C(=O)/C=C\\c1ccccc1",
        "C1=CC=C(C=C1)/C=C/C(=O)/C=C/C2=CC=CC=C2",                                         // dibenzylideneacetone
        "COC1=C(C=CC(=C1)/C=C/C(=O)CC(=O)/C=C/C2=CC(=C(C=C2)O)OC)O",                         // curcumin
        "c1ccccc1c1ccccc1", "[H]/C(F)=C/F", "F/C(/[H])=C/F", "OC(=O)C[C@H](O)[C@H](O)[C@@H](O)CC(=O)O",
    )

    test("every corpus SMILES writes the same canonical form in fifty other atom orders") {
        val random = Random(20261022)
        for (s in corpus) {
            val c = canon(s)
            repeat(50) {
                val other = shuffled(s, random)
                withClue("$s written as $other") { canon(other) shouldBe c }
            }
            withClue("$s: the canonical form is canonical") { canon(c) shouldBe c }
        }
    }

    /** A random SMILES: chains, branches, rings, double bonds with direction marks, stereocenters. Many do not read; those are skipped. */
    fun randomSmiles(r: Random): String {
        val atoms = listOf("C", "C", "C", "N", "O", "F", "Cl", "S", "[C@H]", "[C@@H]", "[C@]", "[C@@]", "[N+]", "[13C]", "Br")
        val sb = StringBuilder()
        var branches = 0
        val open = ArrayList<Int>()
        var lastAtom = false
        val n = r.nextInt(2, 22)
        for (i in 0 until n) {
            if (lastAtom) {
                when (r.nextInt(10)) {
                    0, 1 -> sb.append('=')
                    2 -> sb.append('/')
                    3 -> sb.append('\\')
                    4 -> sb.append('#').also { if (r.nextBoolean()) sb.setLength(sb.length - 1) }
                }
            }
            sb.append(atoms[r.nextInt(atoms.size)])
            lastAtom = true
            if (r.nextInt(6) == 0 && open.size < 3) { val d = (1..9).first { it !in open }; open.add(d); sb.append(d) }
            else if (open.isNotEmpty() && r.nextInt(4) == 0) { sb.append(open.removeAt(r.nextInt(open.size))) }
            if (r.nextInt(5) == 0 && i < n - 1) { sb.append('('); branches++; lastAtom = false.also { sb.append(atoms[r.nextInt(3)]); } ; lastAtom = true }
            if (branches > 0 && r.nextInt(3) == 0) { sb.append(')'); branches-- }
        }
        repeat(branches) { sb.append(')') }
        for (d in open) sb.append("C").append(d)
        return sb.toString()
    }

    test("random molecules with stereochemistry write one canonical form in any atom order (fuzz)") {
        val r = Random(42)
        var tried = 0
        var checked = 0
        var withStereo = 0
        var withDoubleBondStereo = 0
        while (checked < 400 && tried < 20_000) {
            tried++
            val s = randomSmiles(r)
            val mol = Molecule.parseOrNull(s) ?: continue
            checked++
            if (mol.stereo) withStereo++
            if ('/' in mol.smiles || '\\' in mol.smiles) withDoubleBondStereo++
            repeat(20) {
                val other = shuffled(s, r)
                withClue("$s written as $other") { canon(other) shouldBe mol.smiles }
            }
            withClue("$s: idempotent") { canon(mol.smiles) shouldBe mol.smiles }
        }
        checked shouldBe 400
        withClue("$withStereo with stereo, $withDoubleBondStereo with double-bond stereo, of $checked") {
            (withStereo > 100) shouldBe true
        }
    }

    /** A random polyene: double bonds with marked neighbors, conjugated runs, marked branches, rings with marked bonds, stereocenters. */
    fun randomAlkene(r: Random): String {
        fun mark() = if (r.nextBoolean()) "/" else "\\"
        fun sub() = listOf("C", "F", "Cl", "O", "N", "CC", "[13C]", "C(C)C")[r.nextInt(8)]
        val sb = StringBuilder(sub())
        var ring = false
        repeat(r.nextInt(1, 5)) {
            sb.append(mark()).append("C")
            if (!ring && r.nextInt(5) == 0) { sb.append("1"); ring = true }
            if (r.nextInt(3) == 0) sb.append("(").append(mark()).append(sub()).append(")")
            sb.append("=C")
            if (r.nextInt(3) == 0) sb.append("(").append(mark()).append(sub()).append(")")
            if (r.nextInt(4) == 0) sb.append(mark()).append(if (r.nextBoolean()) "[C@H](F)" else "[C@@H](F)")
            if (r.nextInt(4) == 0) sb.append(mark()).append("C(=O)")                                  // a carbonyl between stereo alkenes
        }
        sb.append(mark()).append(sub())
        if (ring) sb.append("CCC1")
        return sb.toString()
    }

    test("random polyenes: double-bond stereo writes one canonical form in any atom order (fuzz)") {
        val r = Random(7)
        var tried = 0
        var checked = 0
        var marked = 0
        while (checked < 300 && tried < 20_000) {
            tried++
            val s = randomAlkene(r)
            val mol = Molecule.parseOrNull(s) ?: continue
            checked++
            if ('/' in mol.smiles || '\\' in mol.smiles) marked++
            repeat(20) {
                val other = try { shuffled(s, r) } catch (e: SciException) { throw AssertionError("$s: ${e.message}") }
                withClue("$s written as $other") { canon(other) shouldBe mol.smiles }
            }
            withClue("$s: idempotent") { canon(mol.smiles) shouldBe mol.smiles }
        }
        checked shouldBe 300
        withClue("$marked of $checked with double-bond stereo") { (marked > 200) shouldBe true }
    }

    test("mirror images: meso compounds equal theirs, chiral ones do not") {
        fun mirror(s: String) = s.replace("@@", "\u0000").replace("@", "@@").replace("\u0000", "@")
        val a = "OC(=O)[C@H](O)[C@H](O)C(=O)O"
        val b = "OC(=O)[C@H](O)[C@@H](O)C(=O)O"
        // Tartaric acid: one of these is meso (its own mirror image), the other chiral.
        listOf(a, b).count { canon(it) == canon(mirror(it)) } shouldBe 1
        canon(a) shouldNotBe canon(b)
        // 1,4-Dimethylcyclohexane: cis and trans, each achiral and different from the other.
        val cis = "C[C@H]1CC[C@@H](C)CC1"
        val trans = "C[C@H]1CC[C@H](C)CC1"
        canon(cis) shouldBe canon(mirror(cis))
        canon(trans) shouldBe canon(mirror(trans))
        canon(cis) shouldNotBe canon(trans)
        // Bromochlorofluoromethane is chiral.
        canon("F[C@H](Cl)Br") shouldNotBe canon(mirror("F[C@H](Cl)Br"))
    }

    test("equal molecules, however written") {
        canon("OCC") shouldBe canon("CCO")
        canon("C(O)C") shouldBe canon("CCO")
        canon("[H]C([H])([H])[H]") shouldBe "C"
        canon("[CH4]") shouldBe "C"
        canon("c1ccccc1") shouldBe canon("C1=CC=CC=C1".let { "c1ccccc1" })
        canon("C1CCCCC1") shouldBe canon("C1CCCCC1")
        canon("C%10CCCCC%10") shouldBe canon("C1CCCCC1")
        canon("[Na+].[Cl-]") shouldBe canon("[Cl-].[Na+]")
        canon("N[C@@H](C)C(=O)O") shouldBe canon("C[C@@H](C(=O)O)N")          // L-alanine, two writings
        canon("N[C@@H](C)C(=O)O") shouldBe canon("OC(=O)[C@@H](N)C")
        canon("F/C=C/F") shouldBe canon("F\\C=C\\F")                           // trans
        canon("F/C=C/F") shouldBe canon("C(\\F)=C/F")
        canon("F/C=C\\F") shouldBe canon("F\\C=C/F")                           // cis
        canon("[H][C@](F)(Cl)Br") shouldBe canon("[C@H](F)(Cl)Br")                  // an explicit H in the same place as the implicit one
        canon("F[C@@H](Cl)Br") shouldBe canon("[C@H](F)(Cl)Br")                    // swapping F and H flips the symbol
        canon("F[C@H](F)Cl") shouldBe canon("F[C@@H](F)Cl")                    // not a real stereocenter
        canon("C=1CCCC=1") shouldBe canon("C1=CCCC1")
        canon("c1ccccc1c1ccccc1") shouldBe canon("c1ccccc1-c1ccccc1")                 // an unwritten bond between rings is single
        canon("[H]/C(F)=C/F") shouldBe canon("F/C=C\\F")                              // a marked [H] folds; its side passes to F
        canon("[H]/C=C/F") shouldBe canon("C=CF")                                      // two hydrogens at one end: no stereo
        canon("C[C@H](C)C") shouldBe canon("CC(C)C")                                   // not stereocenters: the marks go
        canon("F[C@H](F)Cl") shouldBe canon("FC(F)Cl")
        canon("C[C@H]1CCCCC1") shouldBe canon("CC1CCCCC1")
        canon("C/C=C1CCCCC1") shouldBe canon("CC=C1CCCCC1")
        canon("C/C=C/1\\CCCCC1") shouldBe canon("CC=C1CCCCC1")
        canon("C[C@H](C)C[C@H](O)CC") shouldBe canon("CC(C)C[C@H](O)CC")              // one real center among others
        canon("C[C@@H](C)C[C@H](O)CC") shouldBe canon("CC(C)C[C@H](O)CC")
        canon("F/C=C/1.F1") shouldBe canon("F/C=C/F")                          // a ring bond across a dot, marked at its opening
        canon("F/C=C1.F/1") shouldBe canon("F/C=C\\F")                         // marked at its closing: F then C, so C is up from F
    }

    test("different molecules stay different") {
        canon("N[C@@H](C)C(=O)O") shouldNotBe canon("N[C@H](C)C(=O)O")         // L- and D-alanine
        canon("F/C=C/F") shouldNotBe canon("F/C=C\\F")                         // trans and cis
        canon("FC=CF") shouldNotBe canon("F/C=C/F")                            // unspecified and trans
        canon("CCO") shouldNotBe canon("COC")
        canon("c1ccccc1") shouldNotBe canon("C1=CC=CC=C1")                     // no aromaticity model (by design)
        canon("[13CH4]") shouldNotBe canon("C")
        canon("CO") shouldNotBe canon("C[O-]")
        canon("[CH3:1]O") shouldNotBe canon("CO")
        canon("C[S@](=O)CC") shouldNotBe canon("C[S@@](=O)CC")
        canon("[S@](C)(=O)CC") shouldBe canon("C[S@@](=O)CC")                          // a leading lone pair comes first (Open Babel's reading)
        // A pseudoasymmetric center keeps its mark: its two equivalent neighbors hold stereocenters.
        canon("OC(=O)[C@H](O)[C@@H](O)[C@H](O)C(=O)O") shouldNotBe canon("OC(=O)[C@H](O)[C@H](O)[C@H](O)C(=O)O")
        canon("C\\C=C/C(=O)/C=C\\C") shouldNotBe canon("C/C=C/C(=O)/C=C/C")         // (2Z,5Z) and (2E,5E)
        canon("OC[C@H](O)[C@H](O)[C@H](O)CO") shouldNotBe canon("OC[C@H](O)[C@@H](O)[C@H](O)CO")   // ribitol and xylitol: C3 is pseudoasymmetric
    }

    test("canonical forms are plain: brackets only where needed, no bond symbols that are implied") {
        canon("C") shouldBe "C"
        canon("[Na+].[Cl-]") shouldBe "[Cl-].[Na+]"                                   // components sorted by their strings
        canon("c1cc[nH]c1") shouldBe "c1cc[nH]c1"
        canon("OCC") shouldBe "CCO"
        canon("C-C") shouldBe "CC"
        canon("c1ccccc1-c1ccccc1") shouldContain "-"                           // a single bond between aromatic atoms keeps its '-'
        canon("[NH4+]") shouldBe "[NH4+]"
        canon("[Fe+++]") shouldBe "[Fe+3]"
    }

    test("formula, masses and counts") {
        val ethanol = m("CCO")
        ethanol.formula.hill shouldBe "C2H6O"
        near(ethanol.mass, "46.06", "46.08")
        near(ethanol.monoisotopicMass, "46.0418", "46.0419")
        ethanol.atomCount shouldBe 9
        ethanol.heavyAtomCount shouldBe 3
        ethanol.bondCount shouldBe 2
        ethanol.ringCount shouldBe 0
        m("c1ccc2ccccc2c1").ringCount shouldBe 2
        m("[Na+].[Cl-]").components shouldBe 2
        m("[Na+].[Cl-]").charge shouldBe 0
        m("[NH4+]").formula.hill shouldBe "H4N+"
        // PubChem's exact mass of ammonium: the electron is counted.
        near(m("[NH4+]").monoisotopicMass, "18.03382", "18.03383")
        // A labelled atom counts its isotope.
        near(m("[13CH4]").mass, "17.0353", "17.0354")
        m("[2H]O[2H]").formula.hill shouldBe "H2O"
        m("N[C@@H](C)C(=O)O").stereo shouldBe true
        m("CC(=O)O").stereo shouldBe false
        m("c1ccccc1").aromatic shouldBe true
    }

    context("errors say where and why") {
        listOf(
            "" to "Empty", "C(" to "never closed", "C)" to "no '('", "C1CC" to "never closed", "C==C" to "Two bond symbols",
            "=C" to "nothing before", "C=" to "no atom after", "C[Xx]" to "Unknown element", "Na" to "must be written in brackets",
            "Q" to "Unknown atom", "C*C" to "Wildcard", "C C" to "A space", "[C" to "Expected ']'", "CC1C1" to "Two bonds between the same atoms",
            "C11" to "to itself", "Cc" to "outside a ring", "[C@H2](F)Cl" to "cannot be a stereocenter", "[C@](F)Cl" to "neighbors",
            "[999C]" to "no isotope", "[fe]" to "not an aromatic atom", "C[C@OH1](F)Cl" to "Only tetrahedral", "CC)" to "at 3",
            "C=1CC-1" to "two different bond symbols", "F/C=C/1.F/1" to "disagree",
        ).forEach { (smiles, why) ->
            test("'$smiles'") { shouldThrow<SciException> { m(smiles) }.message shouldContain why }
        }
    }

    test("symmetric and large molecules canonicalize quickly") {
        val peptide = "N" + "CC(=O)N[C@@H](C)C(=O)N".repeat(40).dropLast(1) + "O"              // (Gly-Ala)40, about 400 heavy atoms
        val cases = listOf(
            "C(C(C)(C)C)(C(C)(C)C)(C(C)(C)C)C(C)(C)C",                                       // tetra-tert-butylmethane
            "C1C2CC3CC1CC(C2)C3",                                                             // adamantane
            "Cc1c(C)c(C)c(C)c(C)c1C",                                                         // hexamethylbenzene
            "C12C3C4C1C5C2C3C45",                                                             // cubane
            "C1CC2CCC3CCCC4CCC(C1)C2C34",                                                     // a symmetric tricycle
            peptide,
        )
        for (s in cases) {
            val start = System.nanoTime()
            val c = canon(s)
            val ms = (System.nanoTime() - start) / 1_000_000
            withClue("$s took $ms ms") { (ms < 5_000) shouldBe true }
            canon(shuffled(s, Random(3))) shouldBe c
        }
    }

    test("a stored canonical form reads back, for every corpus molecule") {
        for (s in corpus) {
            val m = m(s)
            withClue(s) {
                MoleculeType.decode(MoleculeType.encode(m), 1) shouldBe m
                MoleculeType.parse(MoleculeType.format(m)) shouldBe m
            }
        }
    }

    test("many identical components and long chains are quick") {
        for (s in listOf(List(80) { "C" }.joinToString("."), List(60) { "[Na+].[Cl-]" }.joinToString("."), List(40) { "O" }.joinToString("."))) {
            val start = System.nanoTime()
            val c = canon(s)
            withClue("${s.take(20)}… took ${(System.nanoTime() - start) / 1_000_000} ms") { ((System.nanoTime() - start) / 1_000_000 < 3_000) shouldBe true }
            canon(shuffled(s, Random(5))) shouldBe c
        }
        val chain = "C".repeat(SmilesReader.MAX_ATOMS)
        m(chain).heavyAtomCount shouldBe SmilesReader.MAX_ATOMS                        // no stack overflow at the limit
        MoleculeType.decode(MoleculeType.encode(m(chain)), 1).heavyAtomCount shouldBe SmilesReader.MAX_ATOMS
    }

    test("a molecule too large is refused rather than slow") {
        shouldThrow<SciException> { m("C".repeat(SmilesReader.MAX_ATOMS + 1)) }.message shouldContain "more than"
    }

    test("a long chain and a large symmetric cage canonicalize") {
        m("C".repeat(2000)).atomCount shouldBe 2000 * 3 + 2
        val cubane = m("C12C3C4C1C5C2C3C45")
        cubane.formula.hill shouldBe "C8H8"
        canon(shuffled("C12C3C4C1C5C2C3C45", Random(1))) shouldBe cubane.smiles
    }
})

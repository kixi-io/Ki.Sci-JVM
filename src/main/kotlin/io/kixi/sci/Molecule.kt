package io.kixi.sci

import io.kixi.uom.MolarMass
import io.kixi.uom.Quantity
import java.math.BigDecimal
import io.kixi.uom.Unit as KiUnit

/**
 * A molecular structure, read from SMILES: its atoms (element, isotope, charge,
 * hydrogens, aromatic flag, atom class), bonds and stereochemistry (tetrahedral
 * centers and double bonds).
 *
 * ## Equality
 *
 * Two molecules are equal when their canonical SMILES ([smiles]) are, so
 * `OCC` and `C(O)C` are one molecule, and so are `F/C=C/F` and `F\C=C\F`. The
 * canonical form keeps the structure as written: aromatic benzene (`c1ccccc1`) and
 * a Kekulé form (`C1=CC=CC=C1`) are different values, because Ki.Sci has no
 * aromaticity model (Ki.Chem.Ext, with CDK, will). Explicit hydrogens that carry
 * nothing of their own (`[H]`) are counted on their neighbors, so `[H]C([H])([H])[H]`
 * is `C`.
 *
 * Stereochemistry that is not there is dropped: a center with two equivalent
 * neighbors (`C[C@H](C)C` is `CC(C)C`) and a double bond with two equivalent
 * substituents at one end. An unwritten bond between aromatic atoms outside a ring
 * is single (`c1ccccc1c1ccccc1` is biphenyl), as RDKit, Open Babel and CDK read it.
 * For a three-connected center written first (`[S@](C)(=O)CC`) the lone pair takes
 * the implicit hydrogen's place, first, as Open Babel reads it (RDKit puts it last).
 *
 * Not read: wildcard atoms (`*`), chirality other than tetrahedral (`@`, `@@`,
 * `@TH1`, `@TH2`), and the cis/trans of cumulenes (`F/C=C=C=C/F` reads as
 * `FC=C=C=CF`).
 *
 * ## Masses
 *
 * [mass] is the molar mass from the standard atomic weights (an atom written with an
 * isotope, `[13C]`, counts that isotope's mass), electrons not counted;
 * [monoisotopicMass] uses each element's main isotope and counts the electrons of
 * the net charge, as a mass spectrometer sees the ion.
 *
 * As a value of the `sci` type bag it is written `.smiles("CCO")`, stored as its
 * canonical SMILES, and indexed in the order of that text.
 */
public class Molecule private constructor(
    /** The canonical SMILES. */
    public val smiles: String,
    private val g: SGraph,
) {
    /** The molecular formula: every atom, hydrogens included, and the net charge. Isotopes are not part of a formula. */
    public val formula: Formula by lazy {
        val counts = HashMap<Element, Int>()
        val h = Elements.byNumber(1)!!
        for (a in g.atoms) {
            val e = Elements.byNumber(a.z)!!
            counts[e] = (counts[e] ?: 0) + 1
            if (a.hCount > 0) counts[h] = (counts[h] ?: 0) + a.hCount
        }
        Formula.of(counts, charge)
    }

    /** The net charge. */
    public val charge: Int get() = g.atoms.sumOf { it.charge }

    /** The molar mass, in Da: standard atomic weights, isotope masses where an atom names one, no electrons. */
    public val mass: Quantity<MolarMass> get() = Quantity(massValue, KiUnit.Da)

    internal val massValue: BigDecimal by lazy {
        val h = Elements.byNumber(1)!!.weightValue
        g.atoms.fold(BigDecimal.ZERO) { sum, a ->
            val atom = if (a.isotope != 0) Isotopes.find(a.z, a.isotope)!!.massValue else Elements.byNumber(a.z)!!.weightValue
            sum + atom + h.multiply(BigDecimal(a.hCount))
        }
    }

    /** The monoisotopic mass, in Da: main isotopes (or the isotope an atom names), less the electrons of a positive charge. */
    public val monoisotopicMass: Quantity<MolarMass> get() = Quantity(monoisotopicValue, KiUnit.Da)

    internal val monoisotopicValue: BigDecimal by lazy {
        val h = Elements.byNumber(1)!!.monoisotopicValue
        g.atoms.fold(BigDecimal.ZERO) { sum, a ->
            val atom = if (a.isotope != 0) Isotopes.find(a.z, a.isotope)!!.massValue else Elements.byNumber(a.z)!!.monoisotopicValue
            sum + atom + h.multiply(BigDecimal(a.hCount))
        } - Formula.ELECTRON_MASS_U.multiply(BigDecimal(charge))
    }

    /** Every atom, hydrogens included. */
    public val atomCount: Int get() = g.atoms.size + g.atoms.sumOf { it.hCount }

    /** The atoms other than hydrogen. */
    public val heavyAtomCount: Int get() = g.atoms.count { it.z != 1 }

    /** The bonds between written atoms (bonds to counted hydrogens are not included). */
    public val bondCount: Int get() = g.bonds.size

    /** The disconnected parts (`[Na+].[Cl-]` has two). */
    public val components: Int by lazy {
        val parent = IntArray(g.n) { it }
        fun find(x: Int): Int { var r = x; while (parent[r] != r) { parent[r] = parent[parent[r]]; r = parent[r] }; return r }
        for (b in g.bonds) { val x = find(b.a); val y = find(b.b); if (x != y) parent[x] = y }
        (0 until g.n).count { find(it) == it }
    }

    /** The number of rings (the cycle rank: bonds − atoms + components). */
    public val ringCount: Int get() = g.bonds.size - g.n + components

    /** True when any atom is written aromatic. */
    public val aromatic: Boolean get() = g.atoms.any { it.aromatic }

    /** True when the molecule has stereochemistry: a tetrahedral center or a stereo double bond. */
    public val stereo: Boolean get() = g.atoms.any { it.chiral != 0 } || g.stereo.isNotEmpty()

    override fun equals(other: Any?): Boolean = other is Molecule && other.smiles == smiles
    override fun hashCode(): Int = smiles.hashCode()
    override fun toString(): String = smiles

    public companion object {

        /** Reads SMILES and puts it in canonical form. Throws [SciException] saying where and why. */
        @JvmStatic
        public fun parse(smiles: String): Molecule {
            val g = SmilesReader.read(smiles)
            return Molecule(SmilesCanon.canonical(g).smiles, g)
        }

        /** [parse], or null. */
        @JvmStatic
        public fun parseOrNull(smiles: String): Molecule? = try { parse(smiles) } catch (e: SciException) { null }

        /**
         * A molecule from a string this library wrote as canonical SMILES (RevDB's
         * stored form): read, but not searched again. A change to the canonical form
         * raises the bag's version, so stored strings stay canonical for their version.
         */
        internal fun ofCanonical(smiles: String): Molecule = Molecule(smiles, SmilesReader.read(smiles))
    }
}

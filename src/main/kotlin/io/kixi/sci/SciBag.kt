package io.kixi.sci

import io.kixi.text.ParseException
import io.kixi.types.BagFunction
import io.kixi.types.Member
import io.kixi.types.TypeBag
import io.kixi.types.TypeBags
import io.kixi.types.ValueType
import io.kixi.uom.Quantity
import java.math.BigDecimal
import kotlin.reflect.KClass
import io.kixi.uom.Unit as KiUnit

/**
 * The `sci` type bag: the models of matter that Physics, Chem and Bio read (Ki type
 * bags design, 5.2).
 *
 * - `Element` `.element(Na)`, `Isotope` `.isotope(13C)`, `Formula` `.formula(C6H12O6)`,
 *   `Molecule` `.smiles("CCO")`, `InChI` `.inchi("InChI=1S/…")`.
 * - Functions `element(x)`, `isotope(x)`, `formula(x)`, `molecule(x)`, `mass(x)`.
 * - Constants (CODATA 2022) whose units Ki.Core has: `speedOfLight`, `electronMass`,
 *   `protonMass`, `neutronMass`. The others wait for Ki.Core's energy and charge
 *   units.
 *
 * The bag is found on the classpath through `META-INF/services`. Its Kotlin API is
 * [Elements], [Isotopes], [Formula], [Molecule] and [InChI]. Other bags add members
 * to its types through their `extensions` (Ki.Chem.Ext adds `inchiKey` to Molecule).
 */
public class SciBag : TypeBag {
    override val name: String get() = NAME
    override val version: Int get() = VERSION
    override val types: List<ValueType<*>> get() = TYPES
    override val functions: List<BagFunction> get() = FUNCTIONS
    override val constants: Map<String, Any> get() = CONSTANTS

    public companion object {
        public const val NAME: String = "sci"
        public const val VERSION: Int = 1

        private val TYPES: List<ValueType<*>> = listOf(ElementType, IsotopeType, FormulaType, MoleculeType, InChIType)

        private val FUNCTIONS: List<BagFunction> = listOf(
            BagFunction("element", 1..1, "Element") { a -> elementOf(a[0]) },
            BagFunction("isotope", 1..2, "Isotope") { a -> if (a.size == 2) isotopeOf(a[0], a[1]) else isotopeOf(a[0]) },
            BagFunction("formula", 1..1, "Formula") { a -> formulaOf(a[0]) },
            BagFunction("molecule", 1..1, "Molecule") { a -> moleculeOf(a[0]) },
            BagFunction("mass", 1..1, "Quantity") { a -> massOf(a[0]) },
        )

        /** CODATA 2022 (Mohr et al., Rev. Mod. Phys. 97, 025002, 2025), as quantities. */
        private val CONSTANTS: Map<String, Any> = mapOf(
            "speedOfLight" to Quantity(BigDecimal("299792458"), KiUnit.mps),
            "electronMass" to Quantity(BigDecimal("9.1093837139E-31"), KiUnit.kg),
            "protonMass" to Quantity(BigDecimal("1.67262192595E-27"), KiUnit.kg),
            "neutronMass" to Quantity(BigDecimal("1.67492750056E-27"), KiUnit.kg),
        )

        private fun name(x: Any?): String = x?.let { it::class.simpleName } ?: "nil"

        /** An element from a symbol, a name, an atomic number, an isotope or an element. */
        @JvmStatic
        public fun elementOf(x: Any?): Element = when (x) {
            is Element -> x
            is Isotope -> x.element
            is Number -> Elements.byNumber(x.toInt()).takeIf { x.toDouble() == x.toInt().toDouble() } ?: throw SciException("No element number $x")
            is String -> Elements.lookup(x) ?: throw SciException("Unknown element '$x'")
            else -> throw SciException("Expected an element symbol, name or number, got ${name(x)}")
        }

        /** An isotope from text (`13C`, `C-13`, `D`) or an isotope. */
        @JvmStatic
        public fun isotopeOf(x: Any?): Isotope = when (x) {
            is Isotope -> x
            is String -> Isotopes.lookup(x) ?: throw SciException("Unknown isotope '$x'")
            else -> throw SciException("Expected an isotope such as 13C, got ${name(x)}")
        }

        /** The isotope [massNumber] of an element. */
        @JvmStatic
        public fun isotopeOf(element: Any?, massNumber: Any?): Isotope {
            val e = elementOf(element)
            val a = (massNumber as? Number)?.toInt() ?: (massNumber as? String)?.trim()?.toIntOrNull() ?: throw SciException("Expected a mass number, got ${name(massNumber)}")
            return Isotopes.find(e.number, a) ?: throw SciException("There is no isotope $a${e.symbol}")
        }

        /** A formula from text, a molecule (its formula), or a formula. */
        @JvmStatic
        public fun formulaOf(x: Any?): Formula = when (x) {
            is Formula -> x
            is Molecule -> x.formula
            is String -> Formula.parse(x)
            else -> throw SciException("Expected a formula, got ${name(x)}")
        }

        /** A molecule from SMILES text, or a molecule. */
        @JvmStatic
        public fun moleculeOf(x: Any?): Molecule = when (x) {
            is Molecule -> x
            is String -> Molecule.parse(x)
            else -> throw SciException("Expected SMILES text, got ${name(x)}")
        }

        /**
         * The molar mass of a formula, a molecule or an element (its standard atomic
         * weight), or an isotope's mass, in Da. Text is read as a formula, never as
         * SMILES: `mass("CCO")` is C2O; write `mass(molecule("CCO"))`.
         */
        @JvmStatic
        public fun massOf(x: Any?): Quantity<*> = when (x) {
            is Element -> x.atomicWeight
            is Isotope -> x.mass
            is Molecule -> x.mass
            is Formula -> x.mass
            is InChI -> x.formula?.mass ?: throw SciException("This InChI has no formula layer")
            is String -> Formula.parse(x).mass
            else -> throw SciException("Expected a formula, molecule, element or isotope, got ${name(x)}")
        }
    }
}

/** Reads a literal's inside as a string or token, turning Ki.Sci's exceptions into Ki.Core's ParseException. */
private inline fun <T> readLiteral(literal: String, read: (String) -> T): T {
    val text = TypeBags.stringOrToken(TypeBags.inside(literal))
    return try { read(text) } catch (e: SciException) { throw ParseException(e.message) }
}

/** A KD string literal for [text]: quoted, with `\` and `"` escaped. */
private fun quoted(text: String): String = "\"" + text.replace("\\", "\\\\").replace("\"", "\\\"") + "\""

private fun u16(x: Int): ByteArray = byteArrayOf((x ushr 8).toByte(), x.toByte())
private fun u16(b: ByteArray, at: Int): Int = ((b[at].toInt() and 0xFF) shl 8) or (b[at + 1].toInt() and 0xFF)

/** `Element`, `.element(Na)` (a symbol, a name or an atomic number). Orders by atomic number. */
public object ElementType : ValueType<Element> {
    override val name: String = "Element"
    override val kotlinClass: KClass<Element> = Element::class
    override val literals: List<String> = listOf("element")
    override fun parse(literal: String): Element = readLiteral(literal) { Elements.lookup(it) ?: throw SciException("Unknown element '$it'") }
    override fun format(value: Element): String = ".element(${value.symbol})"
    override fun encode(value: Element): ByteArray = byteArrayOf(value.number.toByte())
    override fun decode(bytes: ByteArray, version: Int): Element {
        require(bytes.size == 1) { "A stored element is one byte" }
        return Elements.byNumber(bytes[0].toInt() and 0xFF) ?: throw IllegalArgumentException("No element ${bytes[0].toInt() and 0xFF}")
    }
    override val ordered: Boolean = true
    override fun sortKey(value: Element): ByteArray = byteArrayOf(value.number.toByte())

    override val members: List<Member<Element>> = listOf(
        Member("number", null, "Int") { e, _ -> e.number },
        Member("symbol", null, "String") { e, _ -> e.symbol },
        Member("name", null, "String") { e, _ -> e.name },
        Member("atomicWeight", null, "Quantity") { e, _ -> e.atomicWeight },
        Member("weightUncertainty", null, "Quantity") { e, _ -> e.weightUncertainty },
        Member("weightRange", null, "Range") { e, _ -> e.weightRange },
        Member("abridgedWeight", null, "Quantity") { e, _ -> e.abridgedWeight },
        Member("nominal", null, "Bool") { e, _ -> e.nominal },
        Member("group", null, "Int") { e, _ -> e.group },
        Member("period", null, "Int") { e, _ -> e.period },
        Member("block", null, "String") { e, _ -> e.block },
        Member("isotopes", null, "List") { e, _ -> e.isotopes },
        Member("mainIsotope", null, "Isotope") { e, _ -> e.mainIsotope },
        Member("monoisotopicMass", null, "Quantity") { e, _ -> e.monoisotopicMass },
    )

    override val samples: List<Element> by lazy { listOf(1, 2, 6, 11, 26, 43, 92, 118).map { Elements.byNumber(it)!! } }
}

/** `Isotope`, `.isotope(13C)`. Orders by atomic number, then mass number. */
public object IsotopeType : ValueType<Isotope> {
    override val name: String = "Isotope"
    override val kotlinClass: KClass<Isotope> = Isotope::class
    override val literals: List<String> = listOf("isotope")
    override fun parse(literal: String): Isotope = readLiteral(literal) { Isotopes.lookup(it) ?: throw SciException("Unknown isotope '$it'") }
    override fun format(value: Isotope): String = ".isotope(${value.symbol})"
    override fun encode(value: Isotope): ByteArray = byteArrayOf(value.element.number.toByte()) + u16(value.massNumber)
    override fun decode(bytes: ByteArray, version: Int): Isotope {
        require(bytes.size == 3) { "A stored isotope is three bytes" }
        val z = bytes[0].toInt() and 0xFF
        val a = u16(bytes, 1)
        return Isotopes.find(z, a) ?: throw IllegalArgumentException("No isotope $a of element $z")
    }
    override val ordered: Boolean = true
    override fun sortKey(value: Isotope): ByteArray = encode(value)

    override val members: List<Member<Isotope>> = listOf(
        Member("element", null, "Element") { i, _ -> i.element },
        Member("symbol", null, "String") { i, _ -> i.symbol },
        Member("massNumber", null, "Int") { i, _ -> i.massNumber },
        Member("neutrons", null, "Int") { i, _ -> i.neutrons },
        Member("mass", null, "Quantity") { i, _ -> i.mass },
        Member("abundance", null, "Double") { i, _ -> i.abundance },
        Member("estimated", null, "Bool") { i, _ -> i.estimated },
    )

    override val samples: List<Isotope> by lazy {
        listOf("1H", "2H", "3H", "12C", "13C", "14C", "16O", "32P", "235U", "238U", "294Og").map { Isotopes.lookup(it)!! }
    }
}

/** `Formula`, `.formula(C6H12O6)`, written in Hill order with its charge. Orders by that text. */
public object FormulaType : ValueType<Formula> {
    override val name: String = "Formula"
    override val kotlinClass: KClass<Formula> = Formula::class
    override val literals: List<String> = listOf("formula")
    override fun parse(literal: String): Formula = readLiteral(literal) { Formula.parse(it) }
    override fun format(value: Formula): String = ".formula(${value.hill})"
    override fun encode(value: Formula): ByteArray = value.hill.toByteArray(Charsets.UTF_8)
    override fun decode(bytes: ByteArray, version: Int): Formula = Formula.parse(String(bytes, Charsets.UTF_8))
    override val ordered: Boolean = true
    override fun sortKey(value: Formula): ByteArray = encode(value)

    override val members: List<Member<Formula>> = listOf(
        Member("hill", null, "String") { f, _ -> f.hill },
        Member("counts", null, "Map") { f, _ -> f.counts },
        Member("elements", null, "List") { f, _ -> f.elements },
        Member("atomCount", null, "Int") { f, _ -> f.atomCount },
        Member("charge", null, "Int") { f, _ -> f.charge },
        Member("mass", null, "Quantity") { f, _ -> f.mass },
        Member("monoisotopicMass", null, "Quantity") { f, _ -> f.monoisotopicMass },
        Member("nominal", null, "Bool") { f, _ -> f.nominal },
        Member("composition", null, "Map") { f, _ -> f.composition },
        Member("count", 1..1, "Int") { f, a -> f.count(SciBag.elementOf(a[0])) },
        Member("times", 1..1, "Formula") { f, a ->
            val n = a[0] as? Number ?: throw SciException("times needs a whole number")
            if (n.toDouble() != n.toInt().toDouble()) throw SciException("times needs a whole number, not $n")
            f.times(n.toInt())
        },
        Member("plus", 1..1, "Formula") { f, a -> f + SciBag.formulaOf(a[0]) },
    )

    override val samples: List<Formula> by lazy {
        listOf("H2O", "C6H12O6", "NaCl", "CuSO4·5H2O", "NH4+", "SO4^2-", "Fe^3+", "H", "C60", "[Cu(NH3)4]SO4", "C12H22O11", "TcO4-")
            .map { Formula.parse(it) }
    }
}

/** `Molecule`, `.smiles("CCO")`, written as canonical SMILES. Orders by that text. */
public object MoleculeType : ValueType<Molecule> {
    override val name: String = "Molecule"
    override val kotlinClass: KClass<Molecule> = Molecule::class
    override val literals: List<String> = listOf("smiles")
    override fun parse(literal: String): Molecule = readLiteral(literal) { Molecule.parse(it) }
    override fun format(value: Molecule): String = ".smiles(${quoted(value.smiles)})"
    override fun encode(value: Molecule): ByteArray = value.smiles.toByteArray(Charsets.UTF_8)
    override fun decode(bytes: ByteArray, version: Int): Molecule = Molecule.ofCanonical(String(bytes, Charsets.UTF_8))
    override val ordered: Boolean = true
    override fun sortKey(value: Molecule): ByteArray = encode(value)

    override val members: List<Member<Molecule>> = listOf(
        Member("smiles", null, "String") { m, _ -> m.smiles },
        Member("formula", null, "Formula") { m, _ -> m.formula },
        Member("mass", null, "Quantity") { m, _ -> m.mass },
        Member("monoisotopicMass", null, "Quantity") { m, _ -> m.monoisotopicMass },
        Member("charge", null, "Int") { m, _ -> m.charge },
        Member("atomCount", null, "Int") { m, _ -> m.atomCount },
        Member("heavyAtomCount", null, "Int") { m, _ -> m.heavyAtomCount },
        Member("bondCount", null, "Int") { m, _ -> m.bondCount },
        Member("ringCount", null, "Int") { m, _ -> m.ringCount },
        Member("components", null, "Int") { m, _ -> m.components },
        Member("aromatic", null, "Bool") { m, _ -> m.aromatic },
        Member("stereo", null, "Bool") { m, _ -> m.stereo },
    )

    override val samples: List<Molecule> by lazy {
        listOf(
            "C", "O", "CCO", "c1ccccc1", "C1=CC=CC=C1", "CC(=O)O", "[Na+].[Cl-]", "[13CH4]", "[2H]O[2H]",
            "N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O", "F/C=C/F", "F/C=C\\F", "C/C=C/C=C/C", "C[S@](=O)CC",
            "c1ccc2ccccc2c1", "C1CC2CCC1CC2", "O=C1C=CC(=O)C=C1", "[NH4+]", "CC(C)(C)C", "C1CCCCC1",
            "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "c1cc[nH]c1", "[O-][N+](=O)c1ccccc1", "C#N", "[Cu+2]", "CCCCCCCCCCCCCCCCCCCC",
        ).map { Molecule.parse(it) }
    }
}

/** `InChI`, `.inchi("InChI=1S/…")`, kept as its text. Orders by that text. */
public object InChIType : ValueType<InChI> {
    override val name: String = "InChI"
    override val kotlinClass: KClass<InChI> = InChI::class
    override val literals: List<String> = listOf("inchi")
    override fun parse(literal: String): InChI = readLiteral(literal) { InChI.parse(it) }
    override fun format(value: InChI): String = ".inchi(${quoted(value.text)})"
    override fun encode(value: InChI): ByteArray = value.text.toByteArray(Charsets.US_ASCII)
    override fun decode(bytes: ByteArray, version: Int): InChI = InChI.parse(String(bytes, Charsets.US_ASCII))
    override val ordered: Boolean = true
    override fun sortKey(value: InChI): ByteArray = encode(value)

    override val members: List<Member<InChI>> = listOf(
        Member("text", null, "String") { i, _ -> i.text },
        Member("standard", null, "Bool") { i, _ -> i.standard },
        Member("layers", null, "List") { i, _ -> i.layers },
        Member("formula", null, "Formula") { i, _ -> i.formula },
    )

    override val samples: List<InChI> by lazy {
        listOf(
            "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
            "InChI=1S/H2O/h1H2",
            "InChI=1S/ClH.Na/h1H;/q;+1/p-1",
            "InChI=1S/p+1",
            "InChI=1/C6H6/c1-2-4-6-5-3-1/h1-6H",
        ).map { InChI.parse(it) }
    }
}

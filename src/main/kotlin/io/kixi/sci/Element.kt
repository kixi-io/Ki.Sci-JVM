package io.kixi.sci

import io.kixi.Range
import io.kixi.uom.MolarMass
import io.kixi.uom.Quantity
import java.math.BigDecimal
import io.kixi.uom.Unit as KiUnit

/**
 * A chemical element: one of the 118 named elements, a singleton per atomic number
 * (two elements are equal when they are the same element).
 *
 * ## Atomic weights
 *
 * [atomicWeight] is the standard atomic weight, CIAAW 2021 (Prohaska et al., Pure
 * Appl. Chem. 94, 2022), in Da (numerically g/mol), with its [weightUncertainty].
 * For the fourteen elements whose weight IUPAC publishes as an interval (H, Li, B,
 * C, N, O, Mg, Si, S, Cl, Ar, Br, Tl, Pb) it is the abridged value, and
 * [weightRange] is the interval. [abridgedWeight] is IUPAC's five-figure table,
 * the one LabNexus's lab kit used. An element with no stable isotope has no
 * standard weight: its weights are the mass number of its longest-lived isotope,
 * and [nominal] is true.
 *
 * ## Isotopes
 *
 * [isotopes] are every nuclide in the AME2020 mass evaluation; [mainIsotope] is the
 * most abundant natural one, or for an element with none the longest-lived, which
 * is what [monoisotopicMass] uses.
 */
public class Element private constructor(
    /** The atomic number. */
    public val number: Int,
    /** The IUPAC symbol, case-sensitive (`Co` is cobalt, `CO` is carbon monoxide). */
    public val symbol: String,
    /** The English name (IUPAC spelling: Aluminium, Caesium, Sulfur). */
    public val name: String,
    private val weight: BigDecimal,
    private val uncertainty: BigDecimal?,
    private val abridged: BigDecimal,
    private val rangeLow: BigDecimal?,
    private val rangeHigh: BigDecimal?,
    /** True when the element has no standard atomic weight: its weights are a mass number. */
    public val nominal: Boolean,
    private val nominalMassNumber: Int,
) {
    /** The standard atomic weight, in Da. */
    public val atomicWeight: Quantity<MolarMass> get() = da(weight)

    /** The standard uncertainty of [atomicWeight], in Da, or null for an element with no standard weight. */
    public val weightUncertainty: Quantity<MolarMass>? get() = uncertainty?.let { da(it) }

    /** The interval IUPAC gives as the standard atomic weight, or null when it gives one value. */
    public val weightRange: Range<Quantity<MolarMass>>? get() =
        if (rangeLow == null || rangeHigh == null) null else Range.inclusive(da(rangeLow), da(rangeHigh))

    /** IUPAC's abridged (five-figure) atomic weight, in Da. */
    public val abridgedWeight: Quantity<MolarMass> get() = da(abridged)

    internal val weightValue: BigDecimal get() = weight

    /** The period (row) of the periodic table, 1 to 7. */
    public val period: Int get() = PERIOD_STARTS.indexOfLast { number >= it } + 1

    /**
     * The group, 1 to 18, or null for the lanthanoids (57 to 71) and actinoids (89 to
     * 103), which IUPAC's table sets apart rather than in group 3.
     */
    public val group: Int? get() = groupOf(number)

    /** The block: `s`, `p`, `d` or `f` (the lanthanoids and actinoids are `f`; helium is `s`). */
    public val block: String get() = when {
        number in 57..71 || number in 89..103 -> "f"
        number == 2 -> "s"
        else -> when (groupOf(number)) { 1, 2 -> "s"; in 3..12 -> "d"; else -> "p" }
    }

    /** Every isotope of the element with a mass in AME2020, by mass number. */
    public val isotopes: List<Isotope> get() = Isotopes.of(this)

    /**
     * The most abundant natural isotope, or for an element with none the longest-lived
     * one when AME2020 has its mass; null for the few superheavy elements it lacks.
     */
    public val mainIsotope: Isotope? get() =
        Isotopes.of(this).filter { it.abundance != null }.maxByOrNull { it.abundance!! }
            ?: Isotopes.find(number, nominalMassNumber)

    /**
     * The mass of [mainIsotope], in Da: what a monoisotopic mass uses. For an element
     * AME2020 has no mass for, the mass number of its longest-lived isotope.
     */
    public val monoisotopicMass: Quantity<MolarMass> get() = da(monoisotopicValue)

    internal val monoisotopicValue: BigDecimal by lazy { mainIsotope?.massValue ?: BigDecimal(nominalMassNumber) }

    override fun toString(): String = symbol

    public companion object {
        private val PERIOD_STARTS = intArrayOf(1, 3, 11, 19, 37, 55, 87)

        private fun da(v: BigDecimal): Quantity<MolarMass> = Quantity(v, KiUnit.Da)

        private fun groupOf(z: Int): Int? = when {
            z == 1 -> 1
            z == 2 -> 18
            z <= 10 -> if (z <= 4) z - 2 else z + 8
            z <= 18 -> if (z <= 12) z - 10 else z
            z <= 36 -> z - 18
            z <= 54 -> z - 36
            z in 57..71 || z in 89..103 -> null
            z <= 56 -> z - 54
            z <= 86 -> z - 68
            z <= 88 -> z - 86
            else -> z - 100
        }

        private fun e(n: Int, s: String, name: String, w: String, u: String, abridged: String, lo: String? = null, hi: String? = null) =
            Element(n, s, name, BigDecimal(w), BigDecimal(u), BigDecimal(abridged), lo?.let(::BigDecimal), hi?.let(::BigDecimal), false, 0)

        private fun n(n: Int, s: String, name: String, massNumber: Int) =
            Element(n, s, name, BigDecimal(massNumber), null, BigDecimal(massNumber), null, null, true, massNumber)

        internal val ALL: List<Element> = listOf(
    e(1, "H", "Hydrogen", "1.0080", "0.0002", "1.008", "1.00784", "1.00811"),
    e(2, "He", "Helium", "4.002602", "0.000002", "4.0026"),
    e(3, "Li", "Lithium", "6.94", "0.06", "6.94", "6.938", "6.997"),
    e(4, "Be", "Beryllium", "9.0121831", "0.0000005", "9.0122"),
    e(5, "B", "Boron", "10.81", "0.02", "10.81", "10.806", "10.821"),
    e(6, "C", "Carbon", "12.011", "0.002", "12.011", "12.0096", "12.0116"),
    e(7, "N", "Nitrogen", "14.007", "0.001", "14.007", "14.00643", "14.00728"),
    e(8, "O", "Oxygen", "15.999", "0.001", "15.999", "15.99903", "15.99977"),
    e(9, "F", "Fluorine", "18.998403162", "0.000000005", "18.998"),
    e(10, "Ne", "Neon", "20.1797", "0.0006", "20.180"),
    e(11, "Na", "Sodium", "22.98976928", "0.00000002", "22.990"),
    e(12, "Mg", "Magnesium", "24.305", "0.002", "24.305", "24.304", "24.307"),
    e(13, "Al", "Aluminium", "26.9815384", "0.0000003", "26.982"),
    e(14, "Si", "Silicon", "28.085", "0.001", "28.085", "28.084", "28.086"),
    e(15, "P", "Phosphorus", "30.973761998", "0.000000005", "30.974"),
    e(16, "S", "Sulfur", "32.06", "0.02", "32.06", "32.059", "32.076"),
    e(17, "Cl", "Chlorine", "35.45", "0.01", "35.45", "35.446", "35.457"),
    e(18, "Ar", "Argon", "39.95", "0.16", "39.95", "39.792", "39.963"),
    e(19, "K", "Potassium", "39.0983", "0.0001", "39.098"),
    e(20, "Ca", "Calcium", "40.078", "0.004", "40.078"),
    e(21, "Sc", "Scandium", "44.955907", "0.000004", "44.956"),
    e(22, "Ti", "Titanium", "47.867", "0.001", "47.867"),
    e(23, "V", "Vanadium", "50.9415", "0.0001", "50.942"),
    e(24, "Cr", "Chromium", "51.9961", "0.0006", "51.996"),
    e(25, "Mn", "Manganese", "54.938043", "0.000002", "54.938"),
    e(26, "Fe", "Iron", "55.845", "0.002", "55.845"),
    e(27, "Co", "Cobalt", "58.933194", "0.000003", "58.933"),
    e(28, "Ni", "Nickel", "58.6934", "0.0004", "58.693"),
    e(29, "Cu", "Copper", "63.546", "0.003", "63.546"),
    e(30, "Zn", "Zinc", "65.38", "0.02", "65.38"),
    e(31, "Ga", "Gallium", "69.723", "0.001", "69.723"),
    e(32, "Ge", "Germanium", "72.630", "0.008", "72.630"),
    e(33, "As", "Arsenic", "74.921595", "0.000006", "74.922"),
    e(34, "Se", "Selenium", "78.971", "0.008", "78.971"),
    e(35, "Br", "Bromine", "79.904", "0.003", "79.904", "79.901", "79.907"),
    e(36, "Kr", "Krypton", "83.798", "0.002", "83.798"),
    e(37, "Rb", "Rubidium", "85.4678", "0.0003", "85.468"),
    e(38, "Sr", "Strontium", "87.62", "0.01", "87.62"),
    e(39, "Y", "Yttrium", "88.905838", "0.000002", "88.906"),
    e(40, "Zr", "Zirconium", "91.224", "0.002", "91.224"),
    e(41, "Nb", "Niobium", "92.90637", "0.00001", "92.906"),
    e(42, "Mo", "Molybdenum", "95.95", "0.01", "95.95"),
    n(43, "Tc", "Technetium", 98),
    e(44, "Ru", "Ruthenium", "101.07", "0.02", "101.07"),
    e(45, "Rh", "Rhodium", "102.90549", "0.00002", "102.91"),
    e(46, "Pd", "Palladium", "106.42", "0.01", "106.42"),
    e(47, "Ag", "Silver", "107.8682", "0.0002", "107.87"),
    e(48, "Cd", "Cadmium", "112.414", "0.004", "112.41"),
    e(49, "In", "Indium", "114.818", "0.001", "114.82"),
    e(50, "Sn", "Tin", "118.710", "0.007", "118.71"),
    e(51, "Sb", "Antimony", "121.760", "0.001", "121.76"),
    e(52, "Te", "Tellurium", "127.60", "0.03", "127.60"),
    e(53, "I", "Iodine", "126.90447", "0.00003", "126.90"),
    e(54, "Xe", "Xenon", "131.293", "0.006", "131.29"),
    e(55, "Cs", "Caesium", "132.90545196", "0.00000006", "132.91"),
    e(56, "Ba", "Barium", "137.327", "0.007", "137.33"),
    e(57, "La", "Lanthanum", "138.90547", "0.00007", "138.91"),
    e(58, "Ce", "Cerium", "140.116", "0.001", "140.12"),
    e(59, "Pr", "Praseodymium", "140.90766", "0.00001", "140.91"),
    e(60, "Nd", "Neodymium", "144.242", "0.003", "144.24"),
    n(61, "Pm", "Promethium", 145),
    e(62, "Sm", "Samarium", "150.36", "0.02", "150.36"),
    e(63, "Eu", "Europium", "151.964", "0.001", "151.96"),
    e(64, "Gd", "Gadolinium", "157.25", "0.03", "157.25"),
    e(65, "Tb", "Terbium", "158.925354", "0.000007", "158.93"),
    e(66, "Dy", "Dysprosium", "162.500", "0.001", "162.50"),
    e(67, "Ho", "Holmium", "164.930329", "0.000005", "164.93"),
    e(68, "Er", "Erbium", "167.259", "0.003", "167.26"),
    e(69, "Tm", "Thulium", "168.934219", "0.000005", "168.93"),
    e(70, "Yb", "Ytterbium", "173.045", "0.010", "173.05"),
    e(71, "Lu", "Lutetium", "174.9668", "0.0001", "174.97"),
    e(72, "Hf", "Hafnium", "178.486", "0.006", "178.49"),
    e(73, "Ta", "Tantalum", "180.94788", "0.00002", "180.95"),
    e(74, "W", "Tungsten", "183.84", "0.01", "183.84"),
    e(75, "Re", "Rhenium", "186.207", "0.001", "186.21"),
    e(76, "Os", "Osmium", "190.23", "0.03", "190.23"),
    e(77, "Ir", "Iridium", "192.217", "0.002", "192.22"),
    e(78, "Pt", "Platinum", "195.084", "0.009", "195.08"),
    e(79, "Au", "Gold", "196.966570", "0.000004", "196.97"),
    e(80, "Hg", "Mercury", "200.592", "0.003", "200.59"),
    e(81, "Tl", "Thallium", "204.38", "0.01", "204.38", "204.382", "204.385"),
    e(82, "Pb", "Lead", "207.2", "1.1", "207.2", "206.14", "207.94"),
    e(83, "Bi", "Bismuth", "208.98040", "0.00001", "208.98"),
    n(84, "Po", "Polonium", 209),
    n(85, "At", "Astatine", 210),
    n(86, "Rn", "Radon", 222),
    n(87, "Fr", "Francium", 223),
    n(88, "Ra", "Radium", 226),
    n(89, "Ac", "Actinium", 227),
    e(90, "Th", "Thorium", "232.0377", "0.0004", "232.04"),
    e(91, "Pa", "Protactinium", "231.03588", "0.00001", "231.04"),
    e(92, "U", "Uranium", "238.02891", "0.00003", "238.03"),
    n(93, "Np", "Neptunium", 237),
    n(94, "Pu", "Plutonium", 244),
    n(95, "Am", "Americium", 243),
    n(96, "Cm", "Curium", 247),
    n(97, "Bk", "Berkelium", 247),
    n(98, "Cf", "Californium", 251),
    n(99, "Es", "Einsteinium", 252),
    n(100, "Fm", "Fermium", 257),
    n(101, "Md", "Mendelevium", 258),
    n(102, "No", "Nobelium", 259),
    n(103, "Lr", "Lawrencium", 266),
    n(104, "Rf", "Rutherfordium", 267),
    n(105, "Db", "Dubnium", 268),
    n(106, "Sg", "Seaborgium", 269),
    n(107, "Bh", "Bohrium", 270),
    n(108, "Hs", "Hassium", 269),
    n(109, "Mt", "Meitnerium", 278),
    n(110, "Ds", "Darmstadtium", 281),
    n(111, "Rg", "Roentgenium", 282),
    n(112, "Cn", "Copernicium", 285),
    n(113, "Nh", "Nihonium", 286),
    n(114, "Fl", "Flerovium", 289),
    n(115, "Mc", "Moscovium", 290),
    n(116, "Lv", "Livermorium", 293),
    n(117, "Ts", "Tennessine", 294),
    n(118, "Og", "Oganesson", 294),
        )
    }
}

/** The periodic table: [all] 118 elements and lookups by symbol, name and number. */
public object Elements {

    /** Every element in atomic-number order. */
    @JvmStatic
    public val all: List<Element> get() = Element.ALL

    private val bySymbol: Map<String, Element> = Element.ALL.associateBy { it.symbol }
    private val byName: Map<String, Element> = Element.ALL.associateBy { it.name.lowercase() } +
            mapOf("aluminum" to Element.ALL[12], "cesium" to Element.ALL[54], "sulphur" to Element.ALL[15])

    /** The element with [symbol] (case-sensitive), or null. */
    @JvmStatic
    public fun bySymbol(symbol: String): Element? = bySymbol[symbol]

    /** The element with English [name] (case-insensitive; the US spellings too), or null. */
    @JvmStatic
    public fun byName(name: String): Element? = byName[name.trim().lowercase()]

    /** The element with atomic [number], or null. */
    @JvmStatic
    public fun byNumber(number: Int): Element? = if (number in 1..118) Element.ALL[number - 1] else null

    /** A symbol (as written, case-sensitive), a name (any case) or an atomic number written as text. */
    @JvmStatic
    public fun lookup(text: String): Element? {
        val t = text.trim()
        return bySymbol[t] ?: byName[t.lowercase()] ?: t.toIntOrNull()?.let { byNumber(it) }
    }
}

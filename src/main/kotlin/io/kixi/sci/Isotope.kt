package io.kixi.sci

import io.kixi.uom.MolarMass
import io.kixi.uom.Quantity
import java.math.BigDecimal
import io.kixi.uom.Unit as KiUnit

/**
 * One nuclide: an [element] and a [massNumber], with its atomic mass from the
 * AME2020 evaluation (Wang et al., Chinese Phys. C 45 030003, 2021) and, for a
 * natural isotope, its representative abundance (CIAAW, Isotopic compositions of
 * the elements 2021). Written `13C`; a singleton per nuclide.
 */
public class Isotope internal constructor(
    /** The element. */
    public val element: Element,
    /** The mass number A (protons plus neutrons). */
    public val massNumber: Int,
    internal val massValue: BigDecimal,
    /** True when AME2020 gives the mass as an estimate rather than a measurement. */
    public val estimated: Boolean,
    /** The representative natural abundance, an amount fraction from 0 to 1, or null for a nuclide not found in nature. */
    public val abundance: Double?,
) {
    /** The atomic mass, in Da (u). */
    public val mass: Quantity<MolarMass> get() = Quantity(massValue, KiUnit.Da)

    /** The neutron count, A − Z. */
    public val neutrons: Int get() = massNumber - element.number

    /** `13C`. */
    public val symbol: String get() = "$massNumber${element.symbol}"

    override fun toString(): String = symbol
}

/** Every nuclide AME2020 has a mass for, loaded once from `isotopes.txt`. */
public object Isotopes {

    private val table: Map<Int, List<Isotope>> by lazy { load() }

    /** The isotope [massNumber] of element [number], or null. */
    @JvmStatic
    public fun find(number: Int, massNumber: Int): Isotope? = table[number]?.firstOrNull { it.massNumber == massNumber }

    /** The isotopes of [element], by mass number. */
    @JvmStatic
    public fun of(element: Element): List<Isotope> = table[element.number] ?: emptyList()

    /**
     * Reads `13C`, `C-13`, `C13`, `carbon-13`, and `D` and `T` for deuterium and
     * tritium. Null when the text names no known nuclide.
     */
    @JvmStatic
    public fun lookup(text: String): Isotope? {
        val t = text.trim()
        when (t) { "D" -> return find(1, 2); "T" -> return find(1, 3) }
        LEADING.matchEntire(t)?.let { m -> return Elements.bySymbol(m.groupValues[2])?.let { find(it.number, m.groupValues[1].toInt()) } }
        TRAILING.matchEntire(t)?.let { m -> return Elements.lookup(m.groupValues[1])?.let { find(it.number, m.groupValues[2].toInt()) } }
        return null
    }

    private val LEADING = Regex("""(\d{1,3})([A-Z][a-z]?)""")
    private val TRAILING = Regex("""([A-Za-z]+)-?(\d{1,3})""")

    private fun load(): Map<Int, List<Isotope>> {
        val stream = Isotopes::class.java.getResourceAsStream("isotopes.txt")
            ?: throw IllegalStateException("Ki.Sci's isotope table (io/kixi/sci/isotopes.txt) is missing from the classpath")
        val out = LinkedHashMap<Int, MutableList<Isotope>>()
        stream.bufferedReader(Charsets.UTF_8).useLines { lines ->
            for (line in lines) {
                if (line.isBlank() || line.startsWith("#")) continue
                val p = line.split(' ')
                val z = p[0].toInt()
                val element = Elements.byNumber(z) ?: continue
                out.getOrPut(z) { ArrayList() }.add(Isotope(element, p[1].toInt(), BigDecimal(p[2]), p[3] == "1",
                    p[4].takeIf { it != "-" }?.toDouble()))
            }
        }
        return out.mapValues { (_, v) -> v.sortedBy { it.massNumber } }
    }
}

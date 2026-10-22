package io.kixi.sci

import io.kixi.uom.MolarMass
import io.kixi.uom.Quantity
import java.math.BigDecimal
import java.math.MathContext
import io.kixi.uom.Unit as KiUnit

/**
 * An empirical or molecular formula: how many atoms of each element, and a net
 * charge. Two formulas are equal when their counts and charges are (`H2O` and `OH2`
 * are one formula; `CO` and `Co` are not). There are no isotopes in a formula; a
 * labelled compound is a [Molecule] (`[13CH4]`).
 *
 * ## Reading
 *
 * - Element symbols are case-sensitive: `Co` is cobalt, `CO` carbon monoxide.
 * - Counts follow what they count: `H2O`, `C6H12O6`; subscript digits (`H₂O`) too.
 * - Groups: `Ca(OH)2`, `Mg3(PO4)2`, `[Cu(NH3)4]SO4`, nested.
 * - Hydrates and adducts join with `·`, `.` or `*`, each part with an optional
 *   leading multiplier: `CuSO4·5H2O`. A leading multiplier on the whole formula
 *   (`2H2O`) multiplies it.
 * - A trailing charge: `NH4+`, `Cl-`, and with a magnitude after a caret, `SO4^2-`,
 *   `Fe^3+`, or in superscripts, `SO₄²⁻`. Written `Fe3+`, the 3 is a count: a
 *   magnitude needs the caret.
 * - Whitespace is ignored.
 *
 * ## Writing
 *
 * In Hill order: carbon, then hydrogen, then the other elements alphabetically;
 * without carbon, every element alphabetically. The charge follows (`H4N+`,
 * `O4S^2-`). This is the literal: `.formula(C6H12O6)`.
 *
 * ## Masses
 *
 * [mass] is the molar mass from the standard atomic weights; electrons are not
 * counted, as a lab weighs an ion's salt. [monoisotopicMass] uses each element's
 * main isotope and does count the electrons a charge adds or removes, as a mass
 * spectrometer sees the ion.
 */
public class Formula private constructor(
    private val atoms: Map<Element, Int>,
    /** The net charge, in units of the elementary charge. */
    public val charge: Int,
) {
    /** The formula in Hill order with its charge: `C6H12O6`, `H4N+`, `O4S^2-`. */
    public val hill: String = buildString {
        for ((e, n) in atoms) { append(e.symbol); if (n != 1) append(n) }
        append(chargeText(charge))
    }

    /** Element symbol to count, in Hill order. */
    public val counts: Map<String, Int> get() = atoms.entries.associateTo(LinkedHashMap()) { it.key.symbol to it.value }

    /** The elements, in Hill order. */
    public val elements: List<Element> get() = atoms.keys.toList()

    /** The count of [element] (0 when it is absent). */
    public fun count(element: Element): Int = atoms[element] ?: 0

    /** Atoms per formula unit. */
    public val atomCount: Int get() = atoms.values.sum()

    /** The molar mass from the standard atomic weights, in Da (numerically g/mol). */
    public val mass: Quantity<MolarMass> get() = Quantity(massValue, KiUnit.Da)

    internal val massValue: BigDecimal by lazy {
        atoms.entries.fold(BigDecimal.ZERO) { a, (e, n) -> a + e.weightValue.multiply(BigDecimal(n)) }
    }

    /** The monoisotopic mass, in Da: each element's main isotope, less the electrons of a positive charge. */
    public val monoisotopicMass: Quantity<MolarMass> get() = Quantity(monoisotopicValue, KiUnit.Da)

    internal val monoisotopicValue: BigDecimal by lazy {
        atoms.entries.fold(BigDecimal.ZERO) { a, (e, n) -> a + e.monoisotopicValue.multiply(BigDecimal(n)) } -
                ELECTRON_MASS_U.multiply(BigDecimal(charge))
    }

    /** True when the mass rests on an element with no standard atomic weight. */
    public val nominal: Boolean get() = atoms.keys.any { it.nominal }

    /** Each element's share of the mass, as a percent, in Hill order. */
    public val composition: Map<String, Double> get() {
        val total = massValue
        return atoms.entries.associateTo(LinkedHashMap()) { (e, n) ->
            e.symbol to e.weightValue.multiply(BigDecimal(n)).multiply(BigDecimal(100)).divide(total, MathContext.DECIMAL64).toDouble()
        }
    }

    /** This formula [n] times, charge included: `H2O.times(2)` is `H4O2`. */
    public fun times(n: Int): Formula {
        if (n < 1) throw SciException("A formula can be multiplied by a whole number from 1, not $n")
        return Formula(atoms.mapValues { Math.multiplyExact(it.value, n) }, Math.multiplyExact(charge, n))
    }

    /** This formula and [other] together (an adduct, a salt), charges added. */
    public operator fun plus(other: Formula): Formula {
        val m = HashMap(atoms)
        for ((e, n) in other.atoms) m[e] = Math.addExact(m[e] ?: 0, n)
        return of(m, Math.addExact(charge, other.charge))
    }

    override fun equals(other: Any?): Boolean = other is Formula && other.hill == hill
    override fun hashCode(): Int = hill.hashCode()
    override fun toString(): String = hill

    public companion object {

        /** The electron's relative atomic mass, CODATA 2022 (5.485 799 090 441(97) × 10⁻⁴). */
        internal val ELECTRON_MASS_U = BigDecimal("0.0005485799090441")

        /** A formula from element counts (none zero or negative) and a charge, put in Hill order. */
        @JvmStatic
        public fun of(counts: Map<Element, Int>, charge: Int = 0): Formula {
            if (counts.isEmpty()) throw SciException("A formula has at least one element")
            for ((e, n) in counts) if (n < 1) throw SciException("The count of ${e.symbol} must be at least 1, not $n")
            return Formula(hillOrder(counts), charge)
        }

        /** Reads a formula (see the class description). Throws [SciException] with the reason. */
        @JvmStatic
        public fun parse(text: String): Formula = FormulaReader(text).read()

        /** [parse], or null. */
        @JvmStatic
        public fun parseOrNull(text: String): Formula? = try { parse(text) } catch (e: SciException) { null }

        internal fun hillOrder(counts: Map<Element, Int>): Map<Element, Int> {
            val c = Elements.bySymbol("C")!!
            val h = Elements.bySymbol("H")!!
            val keys = if (c in counts) listOf(c) + listOfNotNull(h.takeIf { it in counts }) +
                    counts.keys.filter { it != c && it != h }.sortedBy { it.symbol }
            else counts.keys.sortedBy { it.symbol }
            return keys.associateWithTo(LinkedHashMap()) { counts[it]!! }
        }

        internal fun chargeText(charge: Int): String = when {
            charge == 0 -> ""
            charge == 1 -> "+"
            charge == -1 -> "-"
            charge > 0 -> "^$charge+"
            else -> "^${-charge}-"
        }
    }

    /** The formula grammar, over the text with whitespace removed and subscripts mapped. */
    private class FormulaReader(private val original: String) {
        private val s: String
        private val charge: Int

        init {
            val sb = StringBuilder(original.length)
            for (ch in original) {
                if (ch.isWhitespace()) continue
                val sub = SUBSCRIPTS.indexOf(ch)
                sb.append(if (sub >= 0) '0' + sub else ch)
            }
            var t = sb.toString()
            var q = 0
            // Superscript charge: ²⁻, ⁺, ³⁺.
            val sup = Regex("([⁰¹²³⁴⁵⁶⁷⁸⁹]*)([⁺⁻])$").find(t)
            val caret = Regex("""\^(\d*)([+\-])$""").find(t)
            val bare = Regex("""[+\-]$""").find(t)
            when {
                sup != null -> {
                    val digits = sup.groupValues[1].map { '0' + SUPERSCRIPTS.indexOf(it) }.joinToString("")
                    q = (if (digits.isEmpty()) 1 else digits.toIntOrNull() ?: fail("Charge too large")) * (if (sup.groupValues[2] == "⁺") 1 else -1)
                    t = t.substring(0, sup.range.first)
                }
                caret != null -> {
                    val digits = caret.groupValues[1]
                    q = (if (digits.isEmpty()) 1 else digits.toIntOrNull() ?: fail("Charge too large")) * (if (caret.groupValues[2] == "+") 1 else -1)
                    t = t.substring(0, caret.range.first)
                }
                bare != null -> { q = if (bare.value == "+") 1 else -1; t = t.dropLast(1) }
            }
            if (q == 0 && (sup != null || caret != null)) fail("A charge of zero is written without a sign")
            s = t
            charge = q
        }

        private fun fail(why: String): Nothing = throw SciException("$why in formula '$original'")

        fun read(): Formula {
            if (s.isEmpty()) fail("No elements")
            val total = HashMap<Element, Int>()
            for (part in s.split('·', '.', '*', '∙', '•')) {
                if (part.isEmpty()) fail("An empty part")
                var i = 0
                while (i < part.length && part[i].isDigit()) i++
                val multiplier = if (i == 0) 1 else part.substring(0, i).toIntOrNull()?.takeIf { it > 0 } ?: fail("A bad multiplier '${part.substring(0, i)}'")
                val body = part.substring(i)
                if (body.isEmpty()) fail("Nothing after the multiplier $multiplier")
                val p = Parser(body)
                val counts = p.sequence(null)
                if (p.pos != body.length) fail("Unexpected '${body[p.pos]}'")
                for ((e, n) in counts) total[e] = try { Math.addExact(total[e] ?: 0, Math.multiplyExact(n, multiplier)) } catch (x: ArithmeticException) { fail("A count too large") }
            }
            return Formula(hillOrder(total), charge)
        }

        private inner class Parser(private val b: String) {
            var pos = 0

            fun sequence(closer: Char?): HashMap<Element, Int> {
                val counts = HashMap<Element, Int>()
                while (pos < b.length) {
                    val ch = b[pos]
                    when {
                        ch == '(' || ch == '[' -> {
                            pos++
                            val inner = sequence(if (ch == '(') ')' else ']')
                            val n = count()
                            for ((e, c) in inner) counts[e] = try { Math.addExact(counts[e] ?: 0, Math.multiplyExact(c, n)) } catch (x: ArithmeticException) { fail("A count too large") }
                        }
                        ch == closer -> { pos++; return counts }
                        ch == ')' || ch == ']' -> fail("An unmatched '$ch'")
                        ch in 'A'..'Z' -> {
                            val start = pos++
                            if (pos < b.length && b[pos] in 'a'..'z') pos++
                            val symbol = b.substring(start, pos)
                            val e = Elements.bySymbol(symbol) ?: fail("An unknown element '$symbol'")
                            val n = count()
                            counts[e] = try { Math.addExact(counts[e] ?: 0, n) } catch (x: ArithmeticException) { fail("A count too large") }
                        }
                        ch in 'a'..'z' -> fail("'$ch' cannot start an element symbol (symbols are capitalised: Na, Cl)")
                        else -> fail("Unexpected '$ch'")
                    }
                }
                if (closer != null) fail("A missing '$closer'")
                if (counts.isEmpty()) fail("No elements")
                return counts
            }

            private fun count(): Int {
                val start = pos
                while (pos < b.length && b[pos].isDigit()) pos++
                if (pos == start) return 1
                val n = b.substring(start, pos).toIntOrNull() ?: fail("A count too large")
                if (n == 0) fail("A zero count")
                return n
            }
        }

        companion object {
            const val SUBSCRIPTS = "₀₁₂₃₄₅₆₇₈₉"
            const val SUPERSCRIPTS = "⁰¹²³⁴⁵⁶⁷⁸⁹"
        }
    }
}

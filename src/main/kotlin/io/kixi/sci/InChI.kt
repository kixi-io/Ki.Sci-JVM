package io.kixi.sci

/**
 * An IUPAC International Chemical Identifier, `InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3`,
 * kept as its text after a check of its layer structure: the `InChI=` prefix and
 * version, a formula layer in Hill order, and layers that start with a known prefix
 * letter. Two InChIs are equal when their texts are.
 *
 * Ki.Sci does not compute InChIs or their keys: that needs the InChI library, which
 * Ki.Chem.Ext brings (with CDK), adding `inchiKey` and `inchi(molecule)`.
 */
public class InChI private constructor(
    /** The identifier, whole. */
    public val text: String,
) {
    /** True for a standard InChI (`InChI=1S/…`). */
    public val standard: Boolean get() = text.startsWith("InChI=1S/")

    /** The layers after the version, as written: the formula layer first (empty for a bare proton). */
    public val layers: List<String> get() = text.substringAfter('/').split('/')

    /**
     * The formula of the main layer, as written: the neutral formula of the
     * connected structure, components summed (`C2H6O.H2O` is `C2H8O2`). Null for an
     * InChI with no formula layer (a bare proton). A charge or a protonation layer
     * (`/q`, `/p`) is not applied.
     */
    public val formula: Formula? get() {
        val f = layers.first()
        return if (f.isEmpty() || f[0] in PREFIXES) null else formulaOf(f)
    }

    override fun equals(other: Any?): Boolean = other is InChI && other.text == text
    override fun hashCode(): Int = text.hashCode()
    override fun toString(): String = text

    public companion object {
        /** The letters an InChI layer may start with. */
        private const val PREFIXES = "chqpbtmsifor"

        /** Checks [text] as an InChI and keeps it. Throws [SciException] with the reason. */
        @JvmStatic
        public fun parse(text: String): InChI {
            val t = text.trim()
            fun fail(why: String): Nothing = throw SciException("$why in InChI '$t'")
            if (!t.startsWith("InChI=")) fail("An InChI starts with 'InChI='")
            val version = t.removePrefix("InChI=").substringBefore('/')
            if (version != "1" && version != "1S") fail("Version '$version' is not 1 or 1S")
            if (!t.removePrefix("InChI=").contains('/')) fail("No layers")
            if (t.any { it.isWhitespace() || it.code < 0x21 || it.code > 0x7E }) fail("Spaces or characters outside printable ASCII")
            val layers = t.substringAfter('/').split('/')
            val first = layers.first()
            if (first.isEmpty()) fail("An empty formula layer")
            if (first.isNotEmpty() && first[0] !in PREFIXES) {
                try { formulaOf(first) } catch (e: SciException) { fail("The formula layer '$first' is not a formula in Hill order") }
            }
            for ((k, layer) in layers.withIndex()) {
                if (k == 0 && (first.isEmpty() || first[0] !in PREFIXES)) continue
                if (layer.isEmpty()) fail("An empty layer")
                if (layer[0] !in PREFIXES) fail("Layer '/$layer' does not start with a layer letter ($PREFIXES)")
            }
            return InChI(t)
        }

        /** [parse], or null. */
        @JvmStatic
        public fun parseOrNull(text: String): InChI? = try { parse(text) } catch (e: SciException) { null }

        /** The components of a formula layer (`C2H6O.2H2O`), each in Hill order, summed. */
        private fun formulaOf(layer: String): Formula {
            var total: Formula? = null
            for (part in layer.split('.')) {
                val m = Regex("""(\d*)([A-Z][A-Za-z0-9]*)""").matchEntire(part) ?: throw SciException("Bad component '$part'")
                val f = Formula.parse(m.groupValues[2])
                if (f.hill != m.groupValues[2]) throw SciException("Component '$part' is not in Hill order")
                val n = m.groupValues[1].ifEmpty { "1" }.toInt()
                val scaled = if (n == 1) f else f.times(n)
                total = total?.plus(scaled) ?: scaled
            }
            return total!!
        }
    }
}

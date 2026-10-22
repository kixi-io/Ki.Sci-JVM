package io.kixi.sci

/*
 * The SMILES reader (OpenSMILES, without wildcard atoms or the non-tetrahedral
 * chirality classes): text to a molecular graph with its stereochemistry, implicit
 * hydrogens counted and explicit ones folded into their neighbors' counts.
 *
 * The graph is internal: a Molecule is the canonical form (SmilesCanon) of it.
 */

/** An implicit hydrogen in a chiral atom's neighbor order. */
internal const val H_MARK = -1

/** The lone pair of a three-connected chiral atom (a sulfoxide's S), in its neighbor order. */
internal const val LONE_MARK = -2

/** A ring bond not yet closed, in its opening atom's neighbor order (`RING_MARK - number`). */
private const val RING_MARK = -1000

/** Bond orders. AROMATIC is a bond between aromatic atoms. */
internal const val SINGLE = 1
internal const val DOUBLE = 2
internal const val TRIPLE = 3
internal const val QUADRUPLE = 4
internal const val AROMATIC = 5

internal class SAtom(
    val z: Int,
    val aromatic: Boolean,
    val isotope: Int,
    val charge: Int,
    var hCount: Int,
    val bracket: Boolean,
    val atomClass: Int,
    /** 0 none, 1 `@` (anticlockwise), 2 `@@` (clockwise). */
    var chiral: Int,
    /** True when a preceding atom bonds to this one (it is first in [order]). */
    val hasFrom: Boolean,
) {
    /** Neighbors in the order the SMILES met them: atom indices, [H_MARK], [LONE_MARK]. */
    val order: MutableList<Int> = ArrayList(4)

    /** A copy whose neighbor order is renumbered through [map]. */
    fun renumbered(map: IntArray): SAtom = SAtom(z, aromatic, isotope, charge, hCount, bracket, atomClass, chiral, hasFrom).also { c ->
        for (v in order) c.order.add(if (v >= 0) map[v] else v)
    }
}

internal class SBond(var a: Int, var b: Int, var order: Int, val explicit: Boolean = true) {
    /** From direction marks (`/`, `\`): true when [b] is "up" from [a]; null when unmarked. */
    var upB: Boolean? = null
    fun other(x: Int): Int = if (x == a) b else a
}

/**
 * A stereo double bond [d1] = [d2]: the side of each single-bonded substituent,
 * relative to its double-bond atom, in a frame of its own (flipping every side is
 * the same configuration). [side1] is for [d1]'s substituents, [side2] for [d2]'s.
 */
internal class DbStereo(val bond: Int, val d1: Int, val d2: Int, val side1: Map<Int, Boolean>, val side2: Map<Int, Boolean>)

internal class SGraph(val atoms: List<SAtom>, val bonds: List<SBond>, val stereo: List<DbStereo>) {
    val n: Int get() = atoms.size
    /** For each atom, its neighbors and the bonds to them, in bond order. */
    val nbr: Array<IntArray>
    val nbrBond: Array<IntArray>

    init {
        val lists = Array(atoms.size) { ArrayList<Int>() }
        for ((i, b) in bonds.withIndex()) { lists[b.a].add(i); lists[b.b].add(i) }
        nbrBond = Array(atoms.size) { lists[it].toIntArray() }
        nbr = Array(atoms.size) { a -> IntArray(nbrBond[a].size) { k -> bonds[nbrBond[a][k]].other(a) } }
    }

    fun bondBetween(a: Int, b: Int): Int {
        val bs = nbrBond[a]
        for (k in bs.indices) if (nbr[a][k] == b) return bs[k]
        return -1
    }
}

/** The organic subset, written without brackets: symbol to its normal valences. */
internal val ORGANIC_VALENCES: Map<Int, IntArray> = mapOf(
    5 to intArrayOf(3), 6 to intArrayOf(4), 7 to intArrayOf(3, 5), 8 to intArrayOf(2), 15 to intArrayOf(3, 5),
    16 to intArrayOf(2, 4, 6), 9 to intArrayOf(1), 17 to intArrayOf(1), 35 to intArrayOf(1), 53 to intArrayOf(1),
)

/** Atoms that may be written aromatic (lowercase): b, c, n, o, p, s outside brackets; se, as, te inside too. */
private val AROMATIC_ORGANIC = setOf(5, 6, 7, 8, 15, 16)
private val AROMATIC_BRACKET = mapOf("b" to 5, "c" to 6, "n" to 7, "o" to 8, "p" to 15, "s" to 16, "se" to 34, "as" to 33, "te" to 52)

internal fun aromaticSymbolAllowed(z: Int, bracket: Boolean): Boolean =
    z in AROMATIC_ORGANIC || (bracket && z in AROMATIC_BRACKET.values)

/**
 * The hydrogens an organic-subset atom has when written without brackets: the next
 * normal valence above its bond orders, less them; an aromatic atom counts each
 * aromatic bond as one and gives one less (OpenSMILES).
 */
internal fun implicitHydrogens(g: SGraph, a: Int): Int = implicitHydrogens(g.atoms[a], g.nbrBond[a].map { g.bonds[it].order })

internal fun implicitHydrogens(atom: SAtom, orders: List<Int>): Int {
    val valences = ORGANIC_VALENCES[atom.z] ?: return 0
    val sum = orders.sumOf { if (it == AROMATIC) 1 else it }
    val v = valences.firstOrNull { it >= sum } ?: return 0
    val h = v - sum - (if (atom.aromatic) 1 else 0)
    return if (h < 0) 0 else h
}

internal object SmilesReader {

    /** Reads [text] to a graph. Throws [SciException] saying where and why. */
    fun read(text: String): SGraph = Reader(text).read()

    private class RingOpen(val atom: Int, val order: Int, val mark: Int, val placeholder: Int)

    private class Reader(val s: String) {
        val atoms = ArrayList<SAtom>()
        val bonds = ArrayList<SBond>()
        var i = 0

        fun fail(why: String, at: Int = i): Nothing = throw SciException("$why at ${at + 1} in SMILES '$s'")

        fun read(): SGraph {
            if (s.isEmpty()) throw SciException("Empty SMILES")
            var prev = -1
            var pendingOrder = 0
            var pendingMark = 0
            var pendingAt = 0
            val branches = ArrayDeque<Int>()
            val rings = HashMap<Int, RingOpen>()
            var nextRing = 0

            fun bondSymbol(order: Int, mark: Int) {
                if (pendingOrder != 0) fail("Two bond symbols in a row")
                pendingOrder = order; pendingMark = mark; pendingAt = i; i++
            }

            while (i < s.length) {
                when (val c = s[i]) {
                    '(' -> {
                        if (prev < 0) fail("A branch with no atom before it")
                        if (pendingOrder != 0) fail("A bond symbol before '('", pendingAt)
                        branches.addLast(prev); i++
                    }
                    ')' -> {
                        if (branches.isEmpty()) fail("A ')' with no '(' before it")
                        if (pendingOrder != 0) fail("A bond symbol with no atom after it", pendingAt)
                        prev = branches.removeLast(); i++
                    }
                    '.' -> {
                        if (pendingOrder != 0) fail("A bond symbol with no atom after it", pendingAt)
                        prev = -1; i++
                    }
                    '-' -> bondSymbol(SINGLE, 0)
                    '=' -> bondSymbol(DOUBLE, 0)
                    '#' -> bondSymbol(TRIPLE, 0)
                    '$' -> bondSymbol(QUADRUPLE, 0)
                    ':' -> bondSymbol(AROMATIC, 0)
                    '/' -> bondSymbol(SINGLE, 1)
                    '\\' -> bondSymbol(SINGLE, 2)
                    '%', in '0'..'9' -> {
                        if (prev < 0) fail("A ring bond with no atom before it")
                        val at = i
                        val number = if (c == '%') {
                            if (i + 2 >= s.length || !s[i + 1].isDigit() || !s[i + 2].isDigit()) fail("'%' takes two digits (%10 to %99)")
                            s.substring(i + 1, i + 3).toInt().also { i += 3 }
                        } else { i++; c - '0' }
                        val open = rings.remove(number)
                        if (open == null) {
                            val ph = RING_MARK - nextRing++
                            rings[number] = RingOpen(prev, pendingOrder, pendingMark, ph)
                            atoms[prev].order.add(ph)
                        } else {
                            if (open.atom == prev) fail("A ring bond from an atom to itself", at)
                            if (open.order != 0 && pendingOrder != 0 && open.order != pendingOrder) fail("Ring bond $number has two different bond symbols", at)
                            val order = if (open.order != 0) open.order else if (pendingOrder != 0) pendingOrder else defaultOrder(open.atom, prev)
                            val b = addBond(open.atom, prev, order, at, explicit = open.order != 0 || pendingOrder != 0)
                            setMark(b, open.atom, open.mark, at)
                            setMark(b, prev, pendingMark, at)
                            val oo = atoms[open.atom].order
                            oo[oo.indexOf(open.placeholder)] = prev
                            atoms[prev].order.add(open.atom)
                        }
                        pendingOrder = 0; pendingMark = 0
                    }
                    '[' -> prev = addAtom(bracketAtom(prev >= 0), prev, pendingOrder, pendingMark).also { pendingOrder = 0; pendingMark = 0 }
                    '*' -> fail("Wildcard atoms (*) are not supported")
                    else -> {
                        if (!c.isLetter()) fail(if (c.isWhitespace()) "A space" else "Unexpected '$c'")
                        prev = addAtom(organicAtom(prev >= 0), prev, pendingOrder, pendingMark).also { pendingOrder = 0; pendingMark = 0 }
                    }
                }
            }
            if (pendingOrder != 0) fail("A bond symbol with no atom after it", pendingAt)
            if (branches.isNotEmpty()) fail("A '(' that is never closed", s.length - 1)
            if (rings.isNotEmpty()) throw SciException("Ring bond ${rings.keys.sorted().joinToString()} is never closed in SMILES '$s'")
            if (atoms.isEmpty()) throw SciException("No atoms in SMILES '$s'")
            if (atoms.size > MAX_ATOMS) throw SciException("SMILES with more than $MAX_ATOMS atoms are not supported")
            return finish()
        }

        private fun defaultOrder(a: Int, b: Int): Int = if (atoms[a].aromatic && atoms[b].aromatic) AROMATIC else SINGLE

        private fun addAtom(atom: SAtom, prev: Int, pendingOrder: Int, pendingMark: Int): Int {
            val idx = atoms.size
            atoms.add(atom)
            if (prev >= 0) {
                val b = addBond(prev, idx, if (pendingOrder != 0) pendingOrder else defaultOrder(prev, idx), i, explicit = pendingOrder != 0)
                setMark(b, prev, pendingMark, i)
                atoms[prev].order.add(idx)
                atom.order.add(prev)
            } else if (pendingOrder != 0) fail("A bond symbol with nothing before it")
            // An implicit H comes right after the preceding atom (first when there is none).
            if (atom.bracket && atom.hCount == 1 && atom.chiral != 0) atom.order.add(H_MARK)
            return idx
        }

        private val bonded = HashSet<Long>()

        private fun addBond(a: Int, b: Int, order: Int, at: Int, explicit: Boolean): Int {
            if (!bonded.add(minOf(a, b).toLong() shl 32 or maxOf(a, b).toLong())) fail("Two bonds between the same atoms", at)
            bonds.add(SBond(a, b, order, explicit))
            return bonds.size - 1
        }

        /** Records a `/` (1) or `\` (2) written at [from] on bond [b]. */
        private fun setMark(b: Int, from: Int, mark: Int, at: Int) {
            if (mark == 0) return
            val bond = bonds[b]
            val up = if (from == bond.a) mark == 1 else mark == 2
            if (bond.upB != null && bond.upB != up) fail("A ring bond's two direction marks disagree", at)
            bond.upB = up
        }

        private fun organicAtom(hasFrom: Boolean): SAtom {
            val c = s[i]
            if (c.isLowerCase()) {
                val z = when (c) {
                    'b' -> 5; 'c' -> 6; 'n' -> 7; 'o' -> 8; 'p' -> 15; 's' -> 16
                    else -> {
                        // "Na", "Sn": an element outside the organic subset written without brackets.
                        val sym = if (i > 0 && s[i - 1].isUpperCase()) "${s[i - 1]}$c" else null
                        if (sym != null && Elements.bySymbol(sym) != null) fail("$sym must be written in brackets ([$sym])", i - 1)
                        fail("'$c' is not an aromatic atom (b, c, n, o, p, s)")
                    }
                }
                i++
                return SAtom(z, true, 0, 0, 0, false, 0, 0, hasFrom)
            }
            val two = if (i + 1 < s.length) s.substring(i, i + 2) else ""
            val z = when {
                two == "Cl" -> { i += 2; 17 }
                two == "Br" -> { i += 2; 35 }
                else -> {
                    val z1 = when (c) { 'B' -> 5; 'C' -> 6; 'N' -> 7; 'O' -> 8; 'P' -> 15; 'S' -> 16; 'F' -> 9; 'I' -> 53; else -> 0 }
                    if (z1 == 0) {
                        val sym = if (i + 1 < s.length && s[i + 1].isLowerCase()) two else c.toString()
                        fail(if (Elements.bySymbol(sym) != null) "$sym must be written in brackets ([$sym])" else "Unknown atom '$sym'")
                    }
                    i++; z1
                }
            }
            return SAtom(z, false, 0, 0, 0, false, 0, 0, hasFrom)
        }

        private fun bracketAtom(hasFrom: Boolean): SAtom {
            val start = i
            i++ // [
            fun digits(): String { val from = i; while (i < s.length && s[i].isDigit()) i++; return s.substring(from, i) }
            val iso = digits().let { if (it.isEmpty()) 0 else it.toIntOrNull()?.takeIf { v -> v in 1..999 } ?: fail("Isotope $it is out of range", start) }
            if (i >= s.length) fail("An unterminated '['", start)
            var aromatic = false
            val z: Int
            val c = s[i]
            when {
                c == '*' -> fail("Wildcard atoms (*) are not supported")
                c.isUpperCase() -> {
                    val two = if (i + 1 < s.length && s[i + 1].isLowerCase()) s.substring(i, i + 2) else null
                    val e = two?.let { Elements.bySymbol(it) }
                    if (e != null) { z = e.number; i += 2 } else {
                        z = Elements.bySymbol(c.toString())?.number ?: fail("Unknown element '${two ?: c}'")
                        i++
                    }
                }
                c.isLowerCase() -> {
                    val two = if (i + 1 < s.length) s.substring(i, i + 2) else ""
                    val z2 = AROMATIC_BRACKET[two]
                    if (z2 != null && two.length == 2) { z = z2; i += 2 } else {
                        z = AROMATIC_BRACKET[c.toString()] ?: fail("'$c' is not an aromatic atom")
                        i++
                    }
                    aromatic = true
                }
                else -> fail("Expected an element symbol")
            }
            var chiral = 0
            if (i < s.length && s[i] == '@') {
                i++
                chiral = 1
                if (i < s.length && s[i] == '@') { i++; chiral = 2 }
                else if (s.startsWith("TH", i)) {
                    i += 2
                    chiral = when (digits()) { "1" -> 1; "2" -> 2; else -> fail("@TH takes 1 or 2") }
                } else if (i + 1 < s.length && s[i].isUpperCase() && s[i + 1].isUpperCase()) {
                    fail("Only tetrahedral chirality (@, @@) is supported")
                }
            }
            var h = 0
            if (i < s.length && s[i] == 'H') {
                i++
                h = digits().let { if (it.isEmpty()) 1 else it.toIntOrNull()?.takeIf { v -> v <= 9 } ?: fail("Hydrogen count $it is out of range") }
            }
            var charge = 0
            if (i < s.length && (s[i] == '+' || s[i] == '-')) {
                val sign = if (s[i] == '+') 1 else -1
                i++
                val d = digits()
                charge = if (d.isNotEmpty()) sign * (d.toIntOrNull()?.takeIf { it <= 15 } ?: fail("Charge $d is out of range")) else {
                    var k = 1
                    while (i < s.length && s[i] == (if (sign > 0) '+' else '-')) { k++; i++ }
                    sign * k
                }
            }
            var cls = 0
            if (i < s.length && s[i] == ':') {
                i++
                cls = digits().let { if (it.isEmpty()) fail("':' takes an atom class number") else it.toIntOrNull()?.takeIf { v -> v <= 99999 } ?: fail("Atom class $it is out of range") }
            }
            if (i >= s.length || s[i] != ']') fail("Expected ']'")
            i++
            if (aromatic && !aromaticSymbolAllowed(z, true)) fail("${Elements.byNumber(z)!!.symbol} cannot be aromatic", start)
            if (iso != 0 && Isotopes.find(z, iso) == null) fail("There is no isotope $iso${Elements.byNumber(z)!!.symbol}", start)
            return SAtom(z, aromatic, iso, charge, h, true, cls, chiral, hasFrom)
        }

        /** Implicit hydrogens, explicit ones folded, checks, stereo. */
        private fun finish(): SGraph {
            val g0 = SGraph(atoms, bonds, emptyList())
            // An unwritten bond between aromatic atoms outside a ring is single (biphenyl as
            // c1ccccc1c1ccccc1), as RDKit, Open Babel and CDK read it; OpenSMILES would make
            // it aromatic, which outside a ring means nothing.
            val ring = ringBonds(g0)
            for ((b, bond) in bonds.withIndex()) if (bond.order == AROMATIC && !bond.explicit && !ring[b]) bond.order = SINGLE
            for ((a, atom) in atoms.withIndex()) if (!atom.bracket) atom.hCount = implicitHydrogens(g0, a)
            // Stereo double bonds are found before hydrogens fold, so a marked [H] still says which side its neighbor is on.
            val g = foldHydrogens(g0, doubleBondStereo(g0))
            checkAromatic(g)
            for ((a, atom) in g.atoms.withIndex()) checkChiral(g, a, atom)
            return normalizeStereo(g)
        }

        /**
         * Plain explicit hydrogens (`[H]`, no isotope, charge or class, one single bond
         * to an atom that is not hydrogen) become their neighbor's count, so
         * `[H]C([H])([H])[H]` is `C` and `[H]/C(F)=C/F` is `F/C=C\F`. A folded hydrogen
         * leaves the [stereo] double bonds; one left with no substituent at an end (its
         * atom now has two hydrogens) is not a stereo bond.
         */
        private fun foldHydrogens(g: SGraph, stereo: List<DbStereo>): SGraph {
            val drop = BooleanArray(g.n)
            for ((a, atom) in g.atoms.withIndex()) {
                if (atom.z != 1 || atom.isotope != 0 || atom.charge != 0 || atom.atomClass != 0 || atom.chiral != 0 || atom.hCount != 0) continue
                if (g.nbr[a].size != 1) continue
                val bond = g.bonds[g.nbrBond[a][0]]
                val other = g.nbr[a][0]
                if (bond.order != SINGLE || g.atoms[other].z == 1 || drop[other]) continue
                drop[a] = true
                val o = g.atoms[other]
                o.hCount++
                o.order[o.order.indexOf(a)] = H_MARK
            }
            if (drop.none { it }) return SGraph(g.atoms, g.bonds, stereo)
            val map = IntArray(g.n) { -1 }
            val kept = ArrayList<SAtom>()
            for (a in 0 until g.n) if (!drop[a]) { map[a] = kept.size; kept.add(g.atoms[a]) }
            for (atom in kept) for (k in atom.order.indices) { val v = atom.order[k]; if (v >= 0) atom.order[k] = map[v] }
            val bondMap = IntArray(g.bonds.size) { -1 }
            val keptBonds = ArrayList<SBond>()
            for ((i, b) in g.bonds.withIndex()) {
                if (drop[b.a] || drop[b.b]) continue
                bondMap[i] = keptBonds.size
                b.a = map[b.a]; b.b = map[b.b]
                keptBonds.add(b)
            }
            fun remap(side: Map<Int, Boolean>) = side.filterKeys { !drop[it] }.mapKeys { map[it.key] }
            val keptStereo = stereo.mapNotNull { st ->
                val s1 = remap(st.side1)
                val s2 = remap(st.side2)
                if (s1.isEmpty() || s2.isEmpty()) null else DbStereo(bondMap[st.bond], map[st.d1], map[st.d2], s1, s2)
            }
            return SGraph(kept, keptBonds, keptStereo)
        }

        /** Aromatic atoms and bonds must be in rings (a ring bond is any bond that is not a bridge). */
        private fun checkAromatic(g: SGraph) {
            val ring = ringBonds(g)
            for ((b, bond) in g.bonds.withIndex()) {
                if (bond.order == AROMATIC && !ring[b]) throw SciException("An aromatic bond outside a ring in SMILES '$s'")
            }
            for ((a, atom) in g.atoms.withIndex()) {
                if (atom.aromatic && g.nbrBond[a].none { ring[it] }) {
                    throw SciException("An aromatic atom (${Elements.byNumber(atom.z)!!.symbol.lowercase()}) outside a ring in SMILES '$s'")
                }
            }
        }

        /** A tetrahedral center has four neighbors, counting one implicit H or the lone pair of a three-connected atom. */
        private fun checkChiral(g: SGraph, a: Int, atom: SAtom) {
            if (atom.chiral == 0) return
            val explicitNeighbors = g.nbr[a].size
            if (atom.hCount > 1) throw SciException("An atom with ${atom.hCount} hydrogens cannot be a stereocenter in SMILES '$s'")
            val count = explicitNeighbors + atom.hCount
            when (count) {
                4 -> {}
                3 -> if (atom.hCount == 0) atom.order.add(if (atom.hasFrom) 1 else 0, LONE_MARK)
                else throw SciException("Chirality on an atom with three neighbors and a hydrogen is not supported in SMILES '$s'")
                else -> throw SciException("Chirality on an atom with $count neighbors is not supported (tetrahedral centers have 4) in SMILES '$s'")
            }
            if (atom.order.size != 4) throw SciException("Internal error: chiral neighbor order ${atom.order} in SMILES '$s'")
        }

        /**
         * The stereo double bonds: a double bond with a direction-marked single bond at
         * each end. Marks at one end only (a carbonyl between two stereo alkenes, `=CH2`)
         * belong to other double bonds, or to none, and say nothing about this one.
         */
        private fun doubleBondStereo(g: SGraph): List<DbStereo> {
            val out = ArrayList<DbStereo>()
            for ((b, bond) in g.bonds.withIndex()) {
                if (bond.order != DOUBLE) continue
                val e1 = end(g, bond.a, bond.b) ?: continue
                val e2 = end(g, bond.b, bond.a) ?: continue
                out.add(DbStereo(b, bond.a, bond.b, e1.sides(), e2.sides()))
            }
            return out
        }

        /** One end of a double bond: [d]'s single-bonded substituents, and which are marked. */
        private inner class End(val d: Int, val subs: List<Pair<Int, SBond>>) {
            /** The side of each substituent, from the marks; checked, since this end is part of a stereo bond. */
            fun sides(): Map<Int, Boolean> {
                if (subs.size > 2) throw SciException("A double-bond atom with ${subs.size} substituents in SMILES '$s'")
                val marked = subs.filter { it.second.upB != null }
                val side = HashMap<Int, Boolean>()
                for ((x, bond) in marked) side[x] = if (bond.a == d) bond.upB!! else !bond.upB!!
                if (marked.size == 2 && side[subs[0].first] == side[subs[1].first]) {
                    throw SciException("Both substituents of a double-bond atom are marked on the same side in SMILES '$s'")
                }
                if (subs.size == 2 && marked.size == 1) {
                    val m = marked[0].first
                    side[subs.first { it.first != m }.first] = !side[m]!!
                }
                return side
            }
        }

        /** [d]'s end of its double bond to [other], or null when none of its substituents is marked. */
        private fun end(g: SGraph, d: Int, other: Int): End? {
            val subs = ArrayList<Pair<Int, SBond>>()
            for (k in g.nbr[d].indices) {
                val x = g.nbr[d][k]
                val bond = g.bonds[g.nbrBond[d][k]]
                if (x == other || bond.order != SINGLE) continue
                subs.add(x to bond)
            }
            return if (subs.none { it.second.upB != null }) null else End(d, subs)
        }

        /**
         * Drops stereochemistry that is not there: a tetrahedral center with two
         * equivalent neighbors (`C[C@H](C)C`, `F[C@H](F)Cl`) and a double bond with two
         * equivalent substituents at one end (`C/C=C1CCCCC1`), so they equal the same
         * molecule written without it. Two neighbors are equivalent when the refined
         * partition puts them in one cell, they are bonded alike, and neither side of
         * them away from the center holds other stereochemistry (so a pseudoasymmetric
         * center, or cis/trans-1,4-dimethylcyclohexane, keeps its marks). Repeated until
         * nothing changes, as dropping one element can free another.
         */
        private fun normalizeStereo(g: SGraph): SGraph {
            if (g.atoms.none { it.chiral != 0 } && g.stereo.isEmpty()) return g
            val cells = SmilesCanon.refine(g, SmilesCanon.initialCells(g))
            var frames = g.stereo
            var changed = true
            while (changed) {
                changed = false
                fun stereoAtoms(): Set<Int> = (g.atoms.indices.filter { g.atoms[it].chiral != 0 } + frames.flatMap { listOf(it.d1, it.d2) }).toSet()
                for ((a, atom) in g.atoms.withIndex()) {
                    if (atom.chiral == 0) continue
                    if (hasEquivalentPair(g, cells, a, g.nbr[a].toList(), stereoAtoms() - a)) { atom.chiral = 0; changed = true }
                }
                val keep = frames.filter { st ->
                    val others = stereoAtoms() - st.d1 - st.d2
                    !hasEquivalentPair(g, cells, st.d1, st.side1.keys.toList(), others) &&
                            !hasEquivalentPair(g, cells, st.d2, st.side2.keys.toList(), others)
                }
                if (keep.size != frames.size) { frames = keep; changed = true }
            }
            for (atom in g.atoms) if (atom.chiral == 0) atom.order.removeAll { it == H_MARK || it == LONE_MARK }
            return SGraph(g.atoms, g.bonds, frames)
        }

        /** True when two of [center]'s [neighbors] are equivalent and carry no stereochemistry beyond [center]. */
        private fun hasEquivalentPair(g: SGraph, cells: IntArray, center: Int, neighbors: List<Int>, stereo: Set<Int>): Boolean {
            for (i in neighbors.indices) for (j in i + 1 until neighbors.size) {
                val x = neighbors[i]
                val y = neighbors[j]
                if (cells[x] != cells[y]) continue
                if (g.bonds[g.bondBetween(center, x)].order != g.bonds[g.bondBetween(center, y)].order) continue
                if (!reaches(g, x, center, stereo) && !reaches(g, y, center, stereo)) return true
            }
            return false
        }

        /** True when the atoms reachable from [start] without passing [blocked] include one of [targets]. */
        private fun reaches(g: SGraph, start: Int, blocked: Int, targets: Set<Int>): Boolean {
            if (targets.isEmpty()) return false
            val seen = BooleanArray(g.n)
            seen[blocked] = true
            seen[start] = true
            val queue = ArrayDeque<Int>()
            queue.addLast(start)
            while (queue.isNotEmpty()) {
                val v = queue.removeFirst()
                if (v in targets) return true
                for (w in g.nbr[v]) if (!seen[w]) { seen[w] = true; queue.addLast(w) }
            }
            return false
        }
    }

    /** The bonds in rings: every bond that is not a bridge (Tarjan, iteratively). */
    fun ringBonds(g: SGraph): BooleanArray {
        val ring = BooleanArray(g.bonds.size) { true }
        val disc = IntArray(g.n) { -1 }
        val low = IntArray(g.n)
        var time = 0
        for (root in 0 until g.n) {
            if (disc[root] >= 0) continue
            val stack = ArrayDeque<IntArray>()            // [atom, parentBond, next neighbor index]
            disc[root] = time; low[root] = time; time++
            stack.addLast(intArrayOf(root, -1, 0))
            while (stack.isNotEmpty()) {
                val top = stack.last()
                val v = top[0]
                if (top[2] < g.nbr[v].size) {
                    val k = top[2]++
                    val w = g.nbr[v][k]
                    val eb = g.nbrBond[v][k]
                    if (eb == top[1]) continue
                    if (disc[w] < 0) {
                        disc[w] = time; low[w] = time; time++
                        stack.addLast(intArrayOf(w, eb, 0))
                    } else low[v] = minOf(low[v], disc[w])
                } else {
                    stack.removeLast()
                    if (stack.isNotEmpty()) {
                        val p = stack.last()[0]
                        low[p] = minOf(low[p], low[v])
                        if (low[v] > disc[p]) ring[top[1]] = false
                    }
                }
            }
        }
        return ring
    }

    /** The largest molecule read: the canonical search is quadratic in the worst case. */
    const val MAX_ATOMS: Int = 5000
}

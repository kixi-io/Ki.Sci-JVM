package io.kixi.sci

/*
 * Canonical SMILES. Two graphs of the same molecule must write the same string,
 * whatever order their SMILES listed the atoms in.
 *
 * 1. Atoms are partitioned by their invariants (degree, element, aromatic, isotope,
 *    charge, hydrogens, atom class) and the partition is refined until it is
 *    equitable: atoms in one cell see the same cells, by bond order, as each other.
 * 2. While a cell holds more than one atom, the search individualizes each atom of
 *    the first such cell in turn and refines again, down to orders where every atom
 *    has its own rank (the leaves).
 * 3. Each leaf's ranks write a SMILES string; the canonical string is the smallest.
 *    The set of leaves depends only on the molecule, not on the input order, so the
 *    smallest string does too. Stereochemistry is not used to rank; it is written
 *    relative to the leaf's order, so the comparison sees it.
 * 4. Two leaves that write the same string give an automorphism of the molecule
 *    (atom i of one writes where atom i of the other does). A child of a node is
 *    skipped when an automorphism fixing the node's path maps an explored child to
 *    it: its subtree writes the same strings.
 *
 * Changing anything here changes stored sort keys: see the Ki type bags design
 * (Molecule's canonical form is fixed per bag version).
 */
internal object SmilesCanon {

    /** A canonical SMILES and the atoms in the order it writes them. */
    class Result(val smiles: String, val writeOrder: IntArray)

    /** The search gives up beyond this many leaves rather than return a string that may not be canonical. */
    const val MAX_LEAVES: Int = 20_000

    fun canonical(g: SGraph): Result {
        val n = g.n
        val initial = initialCells(g)
        val search = Search(g)
        search.run(initial, IntArray(0))
        return search.best!!.let { Result(it.smiles, it.order) }
    }

    // ------------------------------------------------------------------ partitions

    /** Cell ids by count-rank: an atom's id is the number of atoms in earlier cells. */
    private fun initialCells(g: SGraph): IntArray {
        val inv = Array(g.n) { a ->
            val x = g.atoms[a]
            // Degree first: terminal atoms rank lowest, so a string starts at a chain end (CCO, not C(C)O).
            longArrayOf(g.nbr[a].size.toLong(), x.z.toLong(), if (x.aromatic) 1 else 0, x.isotope.toLong(), x.charge.toLong(),
                x.hCount.toLong(), x.atomClass.toLong())
        }
        return countRank(inv)
    }

    private fun compare(a: LongArray, b: LongArray): Int {
        val m = minOf(a.size, b.size)
        for (k in 0 until m) { val c = a[k].compareTo(b[k]); if (c != 0) return c }
        return a.size.compareTo(b.size)
    }

    private fun countRank(sig: Array<LongArray>): IntArray {
        val n = sig.size
        val idx = (0 until n).sortedWith { x, y -> compare(sig[x], sig[y]) }
        val out = IntArray(n)
        var first = 0
        for (k in 0 until n) {
            if (k > 0 && compare(sig[idx[k]], sig[idx[k - 1]]) != 0) first = k
            out[idx[k]] = first
        }
        return out
    }

    private fun distinct(cells: IntArray): Int = cells.toHashSet().size

    /** Refines [cells] until equitable. Ids stay count-ranks, so an unsplit cell keeps its id. */
    fun refine(g: SGraph, cells: IntArray): IntArray {
        var cur = cells
        var count = distinct(cur)
        val n = g.n
        while (count < n) {
            val sig = Array(n) { a ->
                val s = LongArray(1 + g.nbr[a].size)
                s[0] = cur[a].toLong()
                for (k in g.nbr[a].indices) s[1 + k] = g.bonds[g.nbrBond[a][k]].order.toLong() * n + cur[g.nbr[a][k]]
                java.util.Arrays.sort(s, 1, s.size)
                s
            }
            val next = countRank(sig)
            val c = distinct(next)
            if (c == count) return next
            cur = next
            count = c
        }
        return cur
    }

    // ------------------------------------------------------------------ the search

    private class Leaf(val smiles: String, val order: IntArray)

    private class Search(val g: SGraph) {
        var best: Leaf? = null
        var first: Leaf? = null
        val automorphisms = ArrayList<IntArray>()
        var leaves = 0

        fun run(cells0: IntArray, path: IntArray) {
            val cells = refine(g, cells0)
            val n = g.n
            // The target: the first (lowest id) cell with more than one atom.
            val size = IntArray(n)
            for (c in cells) size[c]++
            var target = -1
            for (c in 0 until n) if (size[c] > 1) { target = c; break }
            if (target < 0) { leaf(cells); return }
            val members = (0 until n).filter { cells[it] == target }
            val explored = ArrayList<Int>()
            for (v in members) {
                if (explored.isNotEmpty() && sameOrbit(v, explored, path)) continue
                explored.add(v)
                val next = cells.copyOf()
                for (a in members) if (a != v) next[a] = target + 1
                run(next, path + v)
            }
        }

        fun leaf(rank: IntArray) {
            if (++leaves > MAX_LEAVES) throw SciException("This molecule is too symmetric to write a canonical SMILES for (over $MAX_LEAVES orders searched)")
            val w = SmilesWriter(g, rank).write()
            val leaf = Leaf(w.smiles, w.writeOrder)
            val f = first
            if (f == null) { first = leaf; best = leaf; return }
            if (leaf.smiles == f.smiles) automorphism(leaf, f)
            val b = best!!
            val c = leaf.smiles.compareTo(b.smiles)
            if (c < 0) best = leaf
            else if (c == 0 && b !== f) automorphism(leaf, b)
        }

        /** The map taking the atom [x] writes at position i to the one [y] writes there. */
        fun automorphism(x: Leaf, y: Leaf) {
            val gamma = IntArray(g.n)
            for (i in x.order.indices) gamma[x.order[i]] = y.order[i]
            if ((0 until g.n).all { gamma[it] == it }) return
            automorphisms.add(gamma)
        }

        /** True when an automorphism fixing every atom of [path] maps an atom of [explored] to [v]. */
        fun sameOrbit(v: Int, explored: List<Int>, path: IntArray): Boolean {
            val usable = automorphisms.filter { gamma -> path.all { gamma[it] == it } }
            if (usable.isEmpty()) return false
            val parent = IntArray(g.n) { it }
            fun find(x: Int): Int { var r = x; while (parent[r] != r) { parent[r] = parent[parent[r]]; r = parent[r] }; return r }
            for (gamma in usable) for (a in 0 until g.n) { val ra = find(a); val rb = find(gamma[a]); if (ra != rb) parent[ra] = rb }
            val rv = find(v)
            return explored.any { find(it) == rv }
        }
    }
}

/**
 * Writes [g] as SMILES with every atom's [rank]: components from their lowest rank,
 * neighbors in rank order, branches before the main chain, ring bonds closed in
 * order and opened with the lowest free number, organic-subset atoms without
 * brackets where their hydrogens are the implicit ones, and stereochemistry
 * relative to the order written.
 */
internal class SmilesWriter(private val g: SGraph, private val rank: IntArray) {

    class Written(val smiles: String, val writeOrder: IntArray)

    private val n = g.n
    private val sortedNbr: Array<IntArray> = Array(n) { a -> g.nbr[a].sortedBy { rank[it] }.toIntArray() }
    private val parent = IntArray(n) { -1 }
    private val parentBond = IntArray(n) { -1 }
    private val emitIndex = IntArray(n) { -1 }
    private val children = Array(n) { ArrayList<Int>() }
    private val ringAt = Array(n) { ArrayList<Int>() }
    private val isRing = BooleanArray(g.bonds.size)
    private val order = ArrayList<Int>(n)
    /** Direction marks to write: bond to "b is up from a". */
    private val marks = HashMap<Int, Boolean>()

    fun write(): Written {
        val starts = (0 until n).sortedBy { rank[it] }
        val sb = StringBuilder()
        for (s in starts) {
            if (emitIndex[s] >= 0) continue
            val compStart = order.size
            dfs(s)
            if (compStart > 0) sb.append('.')
            val compEnd = order.size
            assignMarks(compStart, compEnd)
            emit(s, sb, RingDigits())
        }
        return Written(sb.toString(), order.toIntArray())
    }

    /** Tree edges, ring bonds and the write order (preorder), iteratively. */
    private fun dfs(start: Int) {
        val stack = ArrayDeque<IntArray>()          // [atom, next neighbor index]
        emitIndex[start] = order.size; order.add(start)
        stack.addLast(intArrayOf(start, 0))
        while (stack.isNotEmpty()) {
            val top = stack.last()
            val a = top[0]
            if (top[1] >= sortedNbr[a].size) { stack.removeLast(); continue }
            val v = sortedNbr[a][top[1]++]
            val b = g.bondBetween(a, v)
            if (b == parentBond[a]) continue
            if (emitIndex[v] < 0) {
                parent[v] = a; parentBond[v] = b
                children[a].add(v)
                emitIndex[v] = order.size; order.add(v)
                stack.addLast(intArrayOf(v, 0))
            } else if (!isRing[b] && b != parentBond[v]) {
                isRing[b] = true
                ringAt[a].add(b); ringAt[v].add(b)
            }
        }
    }

    /** Ring bonds at [a] in the order their numbers are written: closings first, then openings, each by partner. */
    private fun ringOrder(a: Int): List<Int> {
        val (closing, opening) = ringAt[a].partition { emitIndex[g.bonds[it].other(a)] < emitIndex[a] }
        return closing.sortedBy { emitIndex[g.bonds[it].other(a)] } + opening.sortedBy { emitIndex[g.bonds[it].other(a)] }
    }

    private class RingDigits {
        val used = BooleanArray(100)
        val digit = HashMap<Int, Int>()
    }

    private fun emit(a: Int, sb: StringBuilder, rings: RingDigits) {
        // Recursion depth is the longest chain; SmilesReader caps the atom count.
        atomText(a, sb)
        val freed = ArrayList<Int>()
        for (b in ringOrder(a)) {
            val other = g.bonds[b].other(a)
            if (emitIndex[other] < emitIndex[a]) {
                val d = rings.digit.remove(b)!!
                sb.append(digitText(d))
                freed.add(d)
            } else {
                sb.append(bondText(b, a))
                val d = (1 until 100).firstOrNull { !rings.used[it] } ?: throw SciException("More than 99 rings open at once")
                rings.used[d] = true
                rings.digit[b] = d
                sb.append(digitText(d))
            }
        }
        for (d in freed) rings.used[d] = false
        val kids = children[a]
        for ((k, c) in kids.withIndex()) {
            val last = k == kids.size - 1
            if (!last) sb.append('(')
            sb.append(bondText(parentBond[c], a))
            emit(c, sb, rings)
            if (!last) sb.append(')')
        }
    }

    private fun digitText(d: Int): String = if (d < 10) d.toString() else "%$d"

    /** The bond symbol for [b], written after atom [from]. */
    private fun bondText(b: Int, from: Int): String {
        val bond = g.bonds[b]
        val x = g.atoms[bond.a]
        val y = g.atoms[bond.b]
        return when (bond.order) {
            SINGLE -> {
                val up = marks[b]
                when {
                    up != null -> { val qUp = if (from == bond.a) up else !up; if (qUp) "/" else "\\" }
                    x.aromatic && y.aromatic -> "-"
                    else -> ""
                }
            }
            DOUBLE -> "="
            TRIPLE -> "#"
            QUADRUPLE -> "$"
            else -> if (x.aromatic && y.aromatic) "" else ":"
        }
    }

    private fun atomText(a: Int, sb: StringBuilder) {
        val atom = g.atoms[a]
        val e = Elements.byNumber(atom.z)!!
        val symbol = if (atom.aromatic) e.symbol.lowercase() else e.symbol
        val chiral = if (atom.chiral == 0) 0 else outputChirality(a, atom)
        val bare = atom.isotope == 0 && atom.charge == 0 && atom.atomClass == 0 && chiral == 0 &&
                ORGANIC_VALENCES.containsKey(atom.z) && (!atom.aromatic || aromaticSymbolAllowed(atom.z, false)) &&
                atom.hCount == implicitHydrogens(g, a)
        if (bare) { sb.append(symbol); return }
        sb.append('[')
        if (atom.isotope != 0) sb.append(atom.isotope)
        sb.append(symbol)
        if (chiral == 1) sb.append('@') else if (chiral == 2) sb.append("@@")
        if (atom.hCount > 0) { sb.append('H'); if (atom.hCount > 1) sb.append(atom.hCount) }
        when {
            atom.charge == 1 -> sb.append('+')
            atom.charge == -1 -> sb.append('-')
            atom.charge > 1 -> sb.append('+').append(atom.charge)
            atom.charge < -1 -> sb.append('-').append(-atom.charge)
        }
        if (atom.atomClass != 0) sb.append(':').append(atom.atomClass)
        sb.append(']')
    }

    /** `@` (1) or `@@` (2) for [a] with its neighbors in the order written here. */
    private fun outputChirality(a: Int, atom: SAtom): Int {
        val out = ArrayList<Int>(4)
        if (parent[a] >= 0) out.add(parent[a])
        if (H_MARK in atom.order) out.add(H_MARK)
        if (LONE_MARK in atom.order) out.add(LONE_MARK)
        for (b in ringOrder(a)) out.add(g.bonds[b].other(a))
        out.addAll(children[a])
        val stored = atom.order
        if (out.size != stored.size || out.toSet() != stored.toSet()) {
            throw SciException("Internal error: chiral neighbors $stored written as $out")
        }
        // Parity of the permutation taking the stored order to the written one.
        val perm = IntArray(out.size) { stored.indexOf(out[it]) }
        var swaps = 0
        val seen = BooleanArray(perm.size)
        for (i in perm.indices) {
            if (seen[i]) continue
            var j = i
            var len = 0
            while (!seen[j]) { seen[j] = true; j = perm[j]; len++ }
            swaps += len - 1
        }
        return if (swaps % 2 == 0) atom.chiral else 3 - atom.chiral
    }

    /**
     * Direction marks for the stereo double bonds of the component written as
     * order[[from], [to]): every single bond at either end of a stereo double bond is
     * marked. A bond between two double bonds (a conjugated chain) carries one mark
     * for both, so their frames are oriented together: the double bonds are taken in
     * write order, and from each one not yet oriented the ones it shares bonds with
     * are reached breadth-first, each oriented by the marks already written. The first
     * double bond of each such group writes its first mark as `/`. Whether a group can
     * be oriented at all does not depend on the write order (only a ring of conjugated
     * stereo double bonds can fail), so an error is the same for every order.
     */
    private fun assignMarks(from: Int, to: Int) {
        val frames = g.stereo.filter { emitIndex[it.d1] in from until to }
            .sortedWith(compareBy({ minOf(emitIndex[it.d1], emitIndex[it.d2]) }, { maxOf(emitIndex[it.d1], emitIndex[it.d2]) }))
        if (frames.isEmpty()) return
        // (bond, double-bond atom, substituent, side) for each frame, bonds in write order.
        val subs = frames.map { st ->
            val list = ArrayList<IntArray>()
            for ((s, side) in st.side1) list.add(intArrayOf(g.bondBetween(st.d1, s), st.d1, s, if (side) 1 else 0))
            for ((s, side) in st.side2) list.add(intArrayOf(g.bondBetween(st.d2, s), st.d2, s, if (side) 1 else 0))
            list.sortWith(compareBy({ minOf(emitIndex[g.bonds[it[0]].a], emitIndex[g.bonds[it[0]].b]) }, { maxOf(emitIndex[g.bonds[it[0]].a], emitIndex[g.bonds[it[0]].b]) }))
            list
        }
        val byBond = HashMap<Int, MutableList<Int>>()
        for ((f, list) in subs.withIndex()) for (e in list) byBond.getOrPut(e[0]) { ArrayList() }.add(f)

        /** The "b up from a" a bond gets in its frame, flipped or not. */
        fun upB(e: IntArray, flip: Boolean): Boolean {
            val sUp = (e[3] == 1) != flip
            return if (g.bonds[e[0]].a == e[1]) sUp else !sUp
        }
        fun conflict(): Nothing = throw SciException("Double-bond stereochemistry that SMILES cannot write consistently (a ring of conjugated stereo double bonds)")

        val done = BooleanArray(frames.size)
        for (root in frames.indices) {
            if (done[root]) continue
            val queue = ArrayDeque<Int>()
            queue.addLast(root)
            done[root] = true
            while (queue.isNotEmpty()) {
                val f = queue.removeFirst()
                var flip: Boolean? = null
                for (e in subs[f]) {
                    val have = marks[e[0]] ?: continue
                    val need = upB(e, false) != have
                    if (flip != null && flip != need) conflict()
                    flip = need
                }
                if (flip == null) {
                    // The first-written mark is "/": the later-written atom up from the earlier one.
                    val e = subs[f][0]
                    val bond = g.bonds[e[0]]
                    val earlierIsA = emitIndex[bond.a] < emitIndex[bond.b]
                    flip = (if (earlierIsA) upB(e, false) else !upB(e, false)) != true
                }
                for (e in subs[f]) {
                    val up = upB(e, flip)
                    val have = marks[e[0]]
                    if (have != null && have != up) conflict()
                    marks[e[0]] = up
                }
                // The frames sharing a bond with this one, in write order.
                val next = subs[f].flatMap { byBond[it[0]]!! }.filter { !done[it] }.distinct().sorted()
                for (x in next) { done[x] = true; queue.addLast(x) }
            }
        }
    }
}

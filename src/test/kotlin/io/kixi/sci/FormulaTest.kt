package io.kixi.sci

import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.doubles.plusOrMinus
import io.kotest.matchers.shouldBe
import io.kotest.matchers.shouldNotBe
import io.kotest.matchers.string.shouldContain
import java.math.BigDecimal

/** The formula grammar, Hill order, masses against textbook values (ported from LabNexus's lab kit tests), and the elements and isotopes. */
class FormulaTest : FunSpec({

    fun mass(text: String): Double = (Formula.parse(text).mass.value as BigDecimal).toDouble()

    context("textbook molar masses") {
        listOf(
            "H2O" to 18.015, "NaCl" to 58.44, "CO2" to 44.009, "C6H12O6" to 180.156, "H2SO4" to 98.079,
            "Ca(OH)2" to 74.093, "Mg3(PO4)2" to 262.86, "CuSO4·5H2O" to 249.68, "(NH4)2SO4" to 132.14, "[Cu(NH3)4]SO4" to 227.73,
        ).forEach { (formula, expected) -> test(formula) { mass(formula) shouldBe (expected plusOrMinus 0.02) } }
    }

    test("full standard weights: NaCl is 58.4398, not the abridged 58.44") {
        Formula.parse("NaCl").mass.toString() shouldBe "58.43976928Da"
    }

    test("hydrate spellings agree") {
        val dot = Formula.parse("CuSO4·5H2O")
        Formula.parse("CuSO4.5H2O") shouldBe dot
        Formula.parse("CuSO4*5H2O") shouldBe dot
        Formula.parse("CuSO4 · 5 H2O") shouldBe dot
        dot.hill shouldBe "CuH10O9S"
    }

    test("Hill order: carbon, hydrogen, then alphabetical; without carbon, all alphabetical") {
        Formula.parse("OH2").hill shouldBe "H2O"
        Formula.parse("HOCH2CH3").hill shouldBe "C2H6O"
        Formula.parse("Ca(OH)2").hill shouldBe "CaH2O2"
        Formula.parse("NaHCO3").hill shouldBe "CHNaO3"
        Formula.parse("2H2O").hill shouldBe "H4O2"
        Formula.parse("ClH").hill shouldBe "ClH"
        Formula.parse("CCl4").counts shouldBe mapOf("C" to 1, "Cl" to 4)
    }

    test("charges: caret, bare sign, superscripts; a count without a caret") {
        Formula.parse("SO4^2-").charge shouldBe -2
        Formula.parse("SO4^2-").hill shouldBe "O4S^2-"
        Formula.parse("NH4+").hill shouldBe "H4N+"
        Formula.parse("Fe^3+").charge shouldBe 3
        Formula.parse("SO₄²⁻") shouldBe Formula.parse("SO4^2-")
        Formula.parse("H₂O") shouldBe Formula.parse("H2O")
        Formula.parse("Fe3+").let { it.count(Elements.bySymbol("Fe")!!) shouldBe 3; it.charge shouldBe 1 }
        Formula.parse("Cl-").charge shouldBe -1
        Formula.parse("NH4+") shouldNotBe Formula.parse("NH4")
        mass("NH4+") shouldBe mass("NH4")                                         // electrons are not weighed
        Formula.parse("NH4+").monoisotopicMass shouldNotBe Formula.parse("NH4").monoisotopicMass
        Formula.parse("SO4 2-") shouldBe Formula.parse("SO4^2-")                    // a charge after a space
        Formula.parse("NH4 +") shouldBe Formula.parse("NH4+")
        Formula.parse("Fe 3+").charge shouldBe 3
    }

    test("symbol case matters") {
        Formula.parse("Co") shouldNotBe Formula.parse("CO")
    }

    test("composition, times and plus") {
        val water = Formula.parse("H2O").composition
        water["H"]!! shouldBe (11.19 plusOrMinus 0.01)
        water["O"]!! shouldBe (88.81 plusOrMinus 0.01)
        Formula.parse("H2O").times(2).hill shouldBe "H4O2"
        (Formula.parse("NaCl") + Formula.parse("H2O")).hill shouldBe "ClH2NaO"
        Formula.parse("SO4^2-").times(2).charge shouldBe -4
        shouldThrow<SciException> { Formula.parse("H2O").times(0) }
    }

    context("errors are readable") {
        listOf("Xx", "nacl", "Ca(OH", "H0", "", "not a compound", "H2O)", "2", "H.", "Fe^0+").forEach { bad ->
            test("'$bad' is refused") { shouldThrow<SciException> { Formula.parse(bad) } }
        }
        test("an unknown element is named") { shouldThrow<SciException> { Formula.parse("NaCk") }.message shouldContain "Ck" }
        test("an overflowing count") { shouldThrow<SciException> { Formula.parse("C2147483647C2") }.message shouldContain "too large" }
        test("groups nested too deep") { shouldThrow<SciException> { Formula.parse("(".repeat(5000) + "H" + ")".repeat(5000)) }.message shouldContain "nested" }
    }

    test("nominal masses") {
        Formula.parse("TcO4-").nominal shouldBe true
        Formula.parse("FeO").nominal shouldBe false
    }

    context("elements") {
        test("the table") {
            Elements.all.size shouldBe 118
            Elements.bySymbol("Fe")!!.name shouldBe "Iron"
            Elements.byName("iron")!!.number shouldBe 26
            Elements.byName("aluminum")!!.symbol shouldBe "Al"
            Elements.lookup("26")!!.symbol shouldBe "Fe"
            Elements.bySymbol("Tc")!!.nominal shouldBe true
            Elements.bySymbol("Na")!!.atomicWeight.toString() shouldBe "22.98976928Da"
            Elements.bySymbol("Na")!!.abridgedWeight.toString() shouldBe "22.99Da"
            Elements.bySymbol("C")!!.weightRange.toString() shouldBe "12.0096Da..12.0116Da"
            Elements.bySymbol("Na")!!.weightRange shouldBe null
        }
        test("groups, periods and blocks") {
            fun e(s: String) = Elements.bySymbol(s)!!
            listOf("H" to 1, "He" to 18, "B" to 13, "Na" to 1, "Ar" to 18, "Fe" to 8, "Hf" to 4, "Og" to 18, "Ra" to 2).forEach { (s, g) -> e(s).group shouldBe g }
            e("La").group shouldBe null
            e("Lu").block shouldBe "f"
            e("He").block shouldBe "s"
            e("Fe").block shouldBe "d"
            e("Cl").block shouldBe "p"
            listOf("H" to 1, "Li" to 2, "Na" to 3, "K" to 4, "Rb" to 5, "Cs" to 6, "Fr" to 7, "Og" to 7).forEach { (s, p) -> e(s).period shouldBe p }
        }
        test("main isotopes and monoisotopic masses") {
            Elements.bySymbol("C")!!.mainIsotope!!.symbol shouldBe "12C"
            Elements.bySymbol("Cl")!!.monoisotopicMass.toString() shouldBe "34.96885269Da"
            Elements.bySymbol("Tc")!!.mainIsotope!!.symbol shouldBe "98Tc"
            Elements.bySymbol("Lr")!!.mainIsotope shouldBe null                    // AME2020 has no mass for 266Lr
            Elements.bySymbol("Lr")!!.monoisotopicMass.toString() shouldBe "266Da"
        }
    }

    context("isotopes") {
        test("spellings") {
            Isotopes.lookup("13C") shouldBe Isotopes.find(6, 13)
            Isotopes.lookup("C-13") shouldBe Isotopes.find(6, 13)
            Isotopes.lookup("C13") shouldBe Isotopes.find(6, 13)
            Isotopes.lookup("carbon-13") shouldBe Isotopes.find(6, 13)
            Isotopes.lookup("D") shouldBe Isotopes.find(1, 2)
            Isotopes.lookup("T")!!.symbol shouldBe "3H"
            Isotopes.lookup("999C") shouldBe null
        }
        test("masses and abundances") {
            val c13 = Isotopes.lookup("13C")!!
            c13.mass.toString() shouldBe "13.00335483534Da"
            c13.abundance!! shouldBe (0.0106 plusOrMinus 1e-9)
            c13.neutrons shouldBe 7
            Isotopes.lookup("14C")!!.abundance shouldBe null
            Elements.bySymbol("U")!!.isotopes.mapNotNull { it.abundance }.sum() shouldBe (1.0 plusOrMinus 1e-9)
        }
    }
})

package io.kixi.sci

import io.kixi.kd.KD
import io.kixi.revdb.RevDB
import io.kixi.types.TypeBags
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.collections.shouldBeEmpty
import io.kotest.matchers.shouldBe
import io.kotest.matchers.string.shouldContain
import io.kotest.matchers.types.shouldBeInstanceOf

/** The `sci` bag as KD, KS, KQ and RevDB see it, and InChI. */
class SciBagTest : FunSpec({

    val bag = TypeBags.byName("sci")!!

    test("the bag is discovered from the classpath, with no problems") {
        bag.shouldBeInstanceOf<SciBag>()
        TypeBags.problems.filter { "sci" in it }.shouldBeEmpty()
        listOf("element", "isotope", "formula", "smiles", "inchi").forEach { TypeBags.literal(it) shouldBe bag.types.first { t -> it in t.literals } }
    }

    test("every sample passes the contract's self-check and the corpus round trip") {
        corpusProblems(bag).shouldBeEmpty()
    }

    test("KD reads the literals, quoted or not") {
        val t = KD.read("""compound name="ethanol" formula=.formula(C2H5OH) structure=.smiles("OCC") bare=.smiles(CCO) el=.element(Na) iso=.isotope(13C) id=.inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")""")
        t["formula"] shouldBe Formula.parse("C2H6O")
        t["structure"] shouldBe Molecule.parse("CCO")
        t["bare"] shouldBe t["structure"]
        t["el"] shouldBe Elements.bySymbol("Na")
        t["iso"] shouldBe Isotopes.lookup("13C")
        t.toString() shouldContain "structure=.smiles(\"CCO\")"
        KD.read("t s=.smiles(\"F/C=C\\\\F\")")["s"] shouldBe Molecule.parse("F/C=C\\F")
        shouldThrow<Exception> { KD.read("t s=.smiles(\"C(\")") }.message!! shouldContain "never closed"
    }

    test("KS: literals, members, functions and constants") {
        ks(".smiles(\"OCC\").formula.hill".let { "say $it" }) shouldBe "C2H6O"
        ks("say .formula(H2O).mass") shouldBe "18.015Da"
        ks("say .element(Fe).name") shouldBe "Iron"
        ks("let m = .smiles(\"CCO\")\nsay m is Molecule\nsay m == .smiles(\"C(O)C\")") shouldBe "true\ntrue"
        ks("use ki.sci.*\nsay mass(\"NaCl\")\nsay element(26).symbol\nsay formula(molecule(\"c1ccccc1\")).hill") shouldBe "58.43976928Da\nFe\nC6H6"
        ks("use ki.sci\nsay sci.speedOfLight") shouldBe "299792458mps"
        ks("say .isotope(13C) < .isotope(14C)") shouldBe "true"
        ks("say .element(Na) < .element(Cl)") shouldBe "true"
    }

    test("KQ: structures and formulas stored, indexed and queried") {
        RevDB.inMemory().use { db ->
            db.writeBlocking {
                createCollection("rev:collection compounds { index structure\n index formula }")
                insert("compounds", "c name=\"ethanol\" structure=.smiles(\"CCO\")")
                insert("compounds", "c name=\"dimethyl ether\" structure=.smiles(\"COC\")")
                insert("compounds", "c name=\"water\" structure=.smiles(\"O\") formula=.formula(H2O)")
            }
            db.kq("from c in compounds where c.structure == .smiles(\"OCC\") select c.name").map { it.values[0] }.toList() shouldBe listOf("ethanol")
            db.kq("from c in compounds where c.structure.formula == .formula(C2H6O) orderby c.name select c.name").map { it.values[0] }.toList() shouldBe
                    listOf("dimethyl ether", "ethanol")
            db.kq("from c in compounds where c.name == \"water\" select c.structure.mass, c.formula.hill").single().values.map { it.toString() } shouldBe
                    listOf("18.015Da", "H2O")
            db.kq("from c in compounds where c.structure.heavyAtomCount > 1 select c.name").size shouldBe 2
            db.kq("from c in compounds select mass(\"NaCl\")").first().values[0].toString() shouldBe "58.43976928Da"
        }
    }

    context("InChI") {
        test("layers and formula") {
            val i = InChI.parse("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
            i.standard shouldBe true
            i.formula!!.hill shouldBe "C2H6O"
            i.layers shouldBe listOf("C2H6O", "c1-2-3", "h3H,2H2,1H3")
            InChI.parse("InChI=1S/C2H6O.H2O/c1-2-3;/h3H,2H2,1H3;1H2").formula!!.hill shouldBe "C2H8O2"
            InChI.parse("InChI=1S/2ClH.Ca/h2*1H;/q;;+2/p-2").formula!!.hill shouldBe "CaCl2H2"
            InChI.parse("InChI=1S/p+1").formula shouldBe null
        }
        test("bad InChIs are refused with the reason") {
            fun bad(s: String, why: String) = shouldThrow<SciException> { InChI.parse(s) }.message shouldContain why
            bad("1S/C2H6O", "starts with 'InChI='")
            bad("InChI=2/C2H6O", "Version")
            bad("InChI=1S", "No layers")
            bad("InChI=1S/C2H6O/x1-2", "layer letter")
            bad("InChI=1S/H6C2O/c1-2-3", "Hill order")
            bad("InChI=1S/C2H6O//c1", "empty layer")
            bad("InChI=1S/C2H6O c1", "Spaces")
            bad("InChI=1S/", "empty formula layer")
        }
    }
})

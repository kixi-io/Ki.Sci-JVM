package io.kixi.sci

import io.kixi.KiException

/** Text Ki.Sci cannot read (a formula, a SMILES string, an InChI), or a value it cannot work with. */
public open class SciException(message: String) : KiException(message)

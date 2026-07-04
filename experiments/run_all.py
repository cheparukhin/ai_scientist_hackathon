import verify_pairs as v
v.COMPOUNDS = {
  # 5-HT2B valvulopathy
  "pergolide":     "CCCN1C[C@@H](C[C@H]2[C@H]1CC3=CNC4=CC=CC2=C34)CSC",
  "fenfluramine":  "CCNC(C)CC1=CC(=CC=C1)C(F)(F)F",
  "norfenfluramine":"CC(CC1=CC(=CC=C1)C(F)(F)F)N",
  "cabergoline":   "CCNC(=O)N(CCCN(C)C)C(=O)[C@@H]1C[C@H]2[C@@H](CC3=CNC4=CC=CC2=C34)N(C1)CC=C",
  "benfluorex":    "CC(CC1=CC(=CC=C1)C(F)(F)F)NCCOC(=O)C2=CC=CC=C2",
  # hERG / QT
  "terfenadine":   "CC(C)(C)C1=CC=C(C=C1)C(CCCN2CCC(CC2)C(C3=CC=CC=C3)(C4=CC=CC=C4)O)O",
  "cisapride":     "CO[C@H]1CN(CC[C@H]1NC(=O)C2=CC(=C(C=C2OC)N)Cl)CCCOC3=CC=C(C=C3)F",
  "astemizole":    "COC1=CC=C(C=C1)CCN2CCC(CC2)NC3=NC4=CC=CC=C4N3CC5=CC=C(C=C5)F",
  "thioridazine":  "CN1CCCCC1CCN2C3=CC=CC=C3SC4=C2C=C(C=C4)SC",
  "sertindole":    "C1CN(CCC1C2=CN(C3=C2C=C(C=C3)Cl)C4=CC=C(C=C4)F)CCN5CCNC5=O",
  "grepafloxacin": "CC1CN(CCN1)C2=C(C(=C3C(=C2)N(C=C(C3=O)C(=O)O)C4CC4)C)F",
}
v.PAIRS = [
  # 5-HT2B valvulopathy (both withdrawn, shared off-target)
  ("pergolide","fenfluramine","5-HT2B agonism -> valve fibrosis; both withdrawn"),
  ("cabergoline","benfluorex","5-HT2B (norfenfluramine metab); both withdrawn"),
  ("pergolide","benfluorex","5-HT2B cross-scaffold"),
  # hERG / QT (both withdrawn, shared pharmacophore)
  ("terfenadine","cisapride","hERG Tyr652/Phe656 cage; both withdrawn"),
  ("astemizole","thioridazine","hERG basic-amine+aromatic; both withdrawn"),
  ("sertindole","grepafloxacin","QT/TdP; sertindole potent hERG, grepaflox weak"),
  ("terfenadine","thioridazine","hERG cross-class"),
  ("astemizole","cisapride","hERG cross-class"),
  # sanity checks (expect HIGHER)
  ("pergolide","cabergoline","SANITY both ergolines -> expect HIGH"),
  ("fenfluramine","norfenfluramine","SANITY parent/metabolite -> expect HIGH"),
]
v.run()

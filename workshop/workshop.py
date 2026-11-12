# /// script
# dependencies = [
#     "marimo",
#     "pandas==3.1.1",
#     "rdkit==2026.9.1",
# ]
# requires-python = ">=3.12"
# ///

import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(PurePath):
    import marimo as mo
    from pathlib import Path

    ROOT: mo.notebook_location()
    DATA_DIR = ROOT / "data"
    return DATA_DIR, mo


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Digital Organic Chemistry Workshop

    In this workshop, we will use **RDKit** to explore how molecular structures can be represented, compared, visualised, and eventually used to make predictions about chemical reactions.

    We will work towards analysing a real high-throughput **Buchwald–Hartwig amination** dataset. But first, we need to get a molecule into Python.

    ## Importing `rdkit`

    The next cell simply imports the necessary packages. Just execute it and don't worry about it.
    """)
    return


@app.cell
def _():
    from rdkit import Chem
    from rdkit.Chem import Draw

    return Chem, Draw


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## From a structure to SMILES

    A **SMILES** string represents a molecular graph using text.

    For example:

    - ethanol: `CCO`
    - acetone: `CC(=O)C`
    - benzene: `c1ccccc1`

    ### Your turn

    Before running the next cells, try writing SMILES strings for:

    1. **propanol**
    2. **cyclohexene**
    3. **bromobenzene**
    4. **_trans_-1,4-dibromo-2-butene** *(challenge)*

    Don't worry if you get them wrong. We will use RDKit to check whether the computer interprets your SMILES as the molecule you intended.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    smiles_input = mo.ui.text(
        value="CCO",
        label="Enter a SMILES string:",
        full_width=True,
    )

    smiles_input
    return (smiles_input,)


@app.cell
def _(Chem, smiles_input):
    smiles = smiles_input.value
    mol = Chem.MolFromSmiles(smiles)
    return mol, smiles


@app.cell(hide_code=True)
def _(Chem, Draw, mo, mol, smiles):
    mo.stop(
        mol is None,
        mo.callout(
            "RDKit could not interpret this SMILES. Check the syntax and try again.",
            kind="warn",
        ),
    )

    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
    image = Draw.MolToImage(mol, size=(300, 200))

    mo.vstack(
        [
            mo.md(f"**You entered:** `{smiles}`"),
            image,
            mo.md(f"**RDKit canonical SMILES:** `{canonical_smiles}`"),
        ],
        align="center",
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## From one molecule to a high-throughput experiment

    We have seen that a molecule can be represented as a **SMILES string** and converted into an RDKit molecular graph.

    Now we will apply exactly the same idea to a real chemical dataset.

    Ahneman *et al.* used high-throughput experimentation to investigate a
    **Buchwald–Hartwig C–N cross-coupling reaction**. Each of the 4,599 experiments contains
    several molecular components together with an experimentally measured reaction yield and some reaction metadata (_e.g._, position on a well-plate).
    """)
    return


@app.cell
def _(DATA_DIR):
    import pandas as pd

    reactions = pd.read_csv(DATA_DIR / "data_table.csv")
    reactions
    return (reactions,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    We will start by looking at the **aryl halides** in the campaign.

    Rather than treating them as rows in a spreadsheet, let's turn their SMILES
    representations back into molecular structures.
    """)
    return


@app.cell
def _(Chem, reactions):
    aryl_halides = (
        reactions[
            ["aryl_halide_number", "aryl_halide", "aryl_halide_smiles"]
        ]
        .drop_duplicates()
        .dropna()
        .sort_values("aryl_halide_number")
        .reset_index(drop=True)
    )

    aryl_halide_mols = [
        Chem.MolFromSmiles(smiles)
        for smiles in aryl_halides["aryl_halide_smiles"]
    ]
 
    aryl_halides
    return aryl_halide_mols, aryl_halides


@app.cell
def _(Chem, Draw, aryl_halide_mols, aryl_halides):
    Draw.MolsToGridImage(
        aryl_halide_mols,
        legends=[
            f"{int(number)}: {Chem.MolToInchiKey(mol)}"
            for number, mol in zip(
                aryl_halides["aryl_halide_number"],
                aryl_halide_mols,
            )
        ],
        molsPerRow=3,
        subImgSize=(280, 190),
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Molecular identifiers

    Each structure above has two identifiers:

    - the **number** used for that aryl halide in the experimental campaign;
    - an **InChIKey** generated from the molecular structure.

    An InChIKey is useful for identifying a molecule, but it does not tell us whether
    two molecules are chemically similar.

    **Look at the structures above.**

    Choose one aryl halide and identify:

    1. another aryl halide that looks **structurally similar**;
    2. another that looks **structurally different**.

    What features did you use to make that judgement?

    Next, we will try to make molecular similarity quantitative.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## From molecular structures to molecular similarity

    Looking at the aryl halides above, we can make intuitive judgements about which
    molecules are structurally similar.

    But how can we make **molecular similarity quantitative**?

    One approach is to convert each molecular graph into a **molecular fingerprint**.

    Here we will use a **Morgan fingerprint**, which describes the local atomic
    environments present in each molecule. We can then compare pairs of fingerprints
    using the **Tanimoto similarity coefficient**.
    """)
    return


@app.cell(hide_code=True)
def _(aryl_halides, mo):
    query_selector = mo.ui.dropdown(
        options={
            f"Aryl halide {int(number)}": i
            for i, number in enumerate(aryl_halides["aryl_halide_number"])
        },
        value="Aryl halide 1",
        label="Choose a query molecule:",
    )
    return (query_selector,)


@app.cell(hide_code=True)
def _(Chem, Draw, aryl_halide_mols, mo, query_selector):
    query_index = query_selector.value
    query_mol = aryl_halide_mols[query_index]

    mo.vstack(
        [
            query_selector,
            Draw.MolToImage(query_mol, size=(300, 220)),
            mo.md(
                f"**InChIKey:** `{Chem.MolToInchiKey(query_mol)}`"
            ),
        ],
        align="center",
    )
    return (query_index,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Morgan fingerprints

    A Morgan fingerprint converts each molecular graph into a fixed-length bit vector
    describing the local atomic environments present in the molecule.

    Here we use:

    - **radius = 2**: environments extending up to two bonds from each atom;
    - **2048 bits**: the length of the fingerprint.
    """)
    return


@app.cell
def _(aryl_halide_mols):
    from rdkit import DataStructs
    from rdkit.Chem import rdFingerprintGenerator

    morgan_generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2,
        fpSize=2048,
    )

    aryl_halide_fps = [
        morgan_generator.GetFingerprint(mol)
        for mol in aryl_halide_mols
    ]
    return DataStructs, aryl_halide_fps


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Once molecules have been converted into binary fingerprints, we can compare
    them using the **Tanimoto similarity coefficient**:

    \[
    T(A,B)=\frac{c}{a+b-c}
    \]

    where:

    - \(a\) is the number of fingerprint features in molecule \(A\);
    - \(b\) is the number of fingerprint features in molecule \(B\);
    - \(c\) is the number of features shared by both molecules.

    Therefore,

    \[
    T(A,B)=
    \frac{\text{shared features}}
    {\text{features present in either molecule}}
    \]

    \(T=1\) means identical fingerprints, while \(T=0\) means that the fingerprints
    have no features in common. No need to implement anything as `rdkit` already has a built in function for such analysis.
    """)
    return


@app.cell
def _(DataStructs, aryl_halide_fps, query_index):
    query_fp = aryl_halide_fps[query_index]

    similarities = [
        DataStructs.TanimotoSimilarity(query_fp, fp)
        for fp in aryl_halide_fps
    ]

    similarity_order = sorted(
        range(len(similarities)),
        key=lambda i: similarities[i],
        reverse=True,
    )

    nearest_indices = similarity_order[:]
    return nearest_indices, similarities


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The aryl halide set **sorted by Tanimoto similarity** to the query molecule is:
    """)
    return


@app.cell(hide_code=True)
def _(Draw, aryl_halide_mols, aryl_halides, nearest_indices, similarities):
    Draw.MolsToGridImage(
        [aryl_halide_mols[i] for i in nearest_indices],
        legends=[
            (
                f"{int(aryl_halides.iloc[i]['aryl_halide_number'])}: "
                f"T = {similarities[i]:.2f}"
            )
            for i in nearest_indices
        ],
        molsPerRow=3,
        subImgSize=(280, 190),
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Molecular descriptors

    Morgan fingerprints encode the **structural environments** present in a molecule.

    But sometimes we want numerical features with a more direct chemical
    interpretation. These are often called **molecular descriptors**.

    For example, in drug discovery, we might describe a molecule using:

    - molecular weight;
    - lipophilicity;
    - polarity;
    - electronic properties.

    Let's calculate a few simple descriptors for the same aryl halides and compare
    them with our fingerprint-based view of similarity.
    """)
    return


@app.cell
def _(aryl_halide_mols, aryl_halides):
    from rdkit.Chem import Descriptors

    aryl_halide_descriptors = aryl_halides[
        ["aryl_halide_number", "aryl_halide_smiles"]
    ].copy()

    aryl_halide_descriptors["MW"] = [
        Descriptors.MolWt(mol)
        for mol in aryl_halide_mols
    ]

    aryl_halide_descriptors["cLogP"] = [
        Descriptors.MolLogP(mol)
        for mol in aryl_halide_mols
    ]

    aryl_halide_descriptors["TPSA"] = [
        Descriptors.TPSA(mol)
        for mol in aryl_halide_mols
    ]

    aryl_halide_descriptors["FractionCSP3"] = [
        Descriptors.FractionCSP3(mol)
        for mol in aryl_halide_mols
    ]

    aryl_halide_descriptors
    return

@app.cell
def _(mo):
    mo.md(r"""
    ### Compare the representations
    
    Look at one of the related **Cl / Br / I** series in the table.
    
    Earlier, the Morgan fingerprints placed these compounds extremely close together.
    
    Now compare their descriptor values.
    
    - Which descriptors change when Cl is replaced by Br or I?
    - Which descriptors remain almost unchanged?
    - Do fingerprints and descriptors therefore contain the same information?
    
    A molecular representation determines **which differences between molecules are
    made visible to a computational model**.
    """)

if __name__ == "__main__":
    app.run()

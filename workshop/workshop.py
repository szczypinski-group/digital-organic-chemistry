# /// script
# dependencies = [
#     "marimo",
#     "pandas==3.1.1",
#     "rdkit==2026.9.1",
# ]
# requires-python = ">=3.14"
# ///

from pathlib import PurePath
import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _():
    import marimo as mo
    from pathlib import Path

    ROOT: PurePath = mo.notebook_location()
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


if __name__ == "__main__":
    app.run()

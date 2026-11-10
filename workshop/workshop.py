# /// script
# dependencies = [
#     "marimo",
#     "rdkit==2026.9.1",
# ]
# requires-python = ">=3.14"
# ///

import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Digital Organic Chemistry Workshop

    In this workshop, we will use **RDKit** to explore how molecular structures can be represented, compared, visualised, and eventually used to make predictions about chemical reactions.

    We will work towards analysing a real high-throughput **Buchwald–Hartwig amination** dataset. But first, we need to get a molecule into Python.



    ### Importing `rdkit`

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
    ## Part 1: from a structure to SMILES

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
    4. **p-bromoanisole** *(challenge)*

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


if __name__ == "__main__":
    app.run()

# /// script
# dependencies = [
#     "altair==6.3.0",
#     "marimo",
#     "numpy==2.5.4",
#     "pandas==3.1.1",
#     "rdkit==2026.9.1",
#     "scikit-learn==1.9.1",
# ]
# requires-python = ">=3.12"
# ///

import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _():
    import marimo as mo
    from pathlib import Path
    import numpy as np

    ROOT: mo.notebook_location()
    # DATA_DIR = ROOT / "data"
    DATA_URL = (
        "https://raw.githubusercontent.com/"
        "szczypinski-group/digital-organic-chemistry/"
        "dev/workshop/data/data_table.csv"
    )
    return DATA_URL, mo, np


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


@app.cell(hide_code=True)
def _(mo):
    mo.image(
        "https://raw.githubusercontent.com/"
        "szczypinski-group/digital-organic-chemistry/"
        "dev/workshop/ahneman-figure.png"
    )
    return


@app.cell(hide_code=True)
def _(DATA_URL):
    import pandas as pd

    reactions = pd.read_csv(DATA_URL)
    reactions
    return pd, reactions


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    We will start by looking at the **aryl halides** in the campaign.

    Rather than treating them as rows in a spreadsheet, let's turn their SMILES
    representations back into molecular structures.
    """)
    return


@app.cell(hide_code=True)
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


@app.cell(hide_code=True)
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
    """)
    return


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


@app.cell(hide_code=True)
def _(mo):
    fp_radius = mo.ui.slider(
        start=0,
        stop=3,
        step=1,
        value=2,
        label="Morgan fingerprint radius",
        show_value=True,
    )

    fp_size = mo.ui.slider(
        start=256,
        stop=2048,
        step=256,
        value=1024,
        label="Fingerprint size",
        show_value=True,
    )

    mo.hstack(
        [fp_radius, fp_size],
        gap=2,
        justify="start",
    )
    return fp_radius, fp_size


@app.cell(hide_code=True)
def _(aryl_halide_mols, fp_radius, fp_size, query_index):
    from rdkit import DataStructs
    from rdkit.Chem import rdFingerprintGenerator

    # FINGERPRINT GENERATION

    morgan_generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=fp_radius.value,
        fpSize=fp_size.value,
    )

    aryl_halide_fps = [
        morgan_generator.GetFingerprint(mol)
        for mol in aryl_halide_mols
    ]

    # TANIMOTO SIMILARITY CALCULATION

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
    return aryl_halide_fps, nearest_indices, similarities


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
    The the **PCA representation of the ECFP space** is:
    """)
    return


@app.cell(hide_code=True)
def _(PCA, alt, aryl_halide_fps, aryl_halides, np, query_index):
    fingerprint_matrix = np.array([
        list(fp)
        for fp in aryl_halide_fps
    ])

    ecfp_pca = PCA(n_components=2)
    ecfp_scores = ecfp_pca.fit_transform(fingerprint_matrix)

    ecfp_pca_data = aryl_halides[
        ["aryl_halide_number", "aryl_halide_smiles"]
    ].copy()

    ecfp_pca_data["PC1"] = ecfp_scores[:, 0]
    ecfp_pca_data["PC2"] = ecfp_scores[:, 1]

    ecfp_pca_data["selected"] = [
        i == query_index
        for i in range(len(ecfp_pca_data))
    ]

    ecfp_variance = ecfp_pca.explained_variance_ratio_ * 100

    ecfp_plot = (
        alt.Chart(ecfp_pca_data)
        .mark_circle(size=45, opacity=0.7)
        .encode(
            x=alt.X(
                "PC1:Q",
                title=f"PC1 ({ecfp_variance[0]:.1f}% variance)",
                scale=alt.Scale(zero=False),
            ),
            y=alt.Y(
                "PC2:Q",
                title=f"PC2 ({ecfp_variance[1]:.1f}% variance)",
                scale=alt.Scale(zero=False),
            ),
            color=alt.condition(
                "datum.selected",
                alt.value("#d62728"),
                alt.value("#777777"),
            ),
            tooltip=[
                alt.Tooltip(
                    "aryl_halide_number:Q",
                    title="Aryl halide",
                    format=".0f",
                ),
                alt.Tooltip(
                    "aryl_halide_smiles:N",
                    title="SMILES",
                ),
            ],
        )
        .properties(
            width=600,
            height=430,
            title="Morgan fingerprint space",
        )
        .interactive()
    )

    ecfp_plot
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Molecular descriptors

    Morgan fingerprints encode the **structural environments** present in a molecule.

    But sometimes we want numerical features with a more direct chemical interpretation. These are often called **molecular descriptors**.

    For the same aryl halides, we will calculate descriptors describing:

    - **size:** molecular weight (MW), heavy-atom count, and Labute accessible surface area;
    - **lipophilicity and polarisability:** cLogP and molecular refractivity (MolMR);
    - **polarity and hydrogen bonding:** topological polar surface area (TPSA), hydrogen-bond acceptors (HBA), and hydrogen-bond donors (HBD);
    - **shape and topology:** rotatable bonds, fraction of sp³ carbon atoms, ring count, and aromatic ring count.

    We will then use PCA to explore how these molecular properties distinguish the aryl halides.
    """)
    return


@app.cell(hide_code=True)
def _(aryl_halide_mols, aryl_halides):
    from rdkit.Chem import Descriptors, Lipinski

    descriptor_data = aryl_halides[
        ["aryl_halide_number", "aryl_halide_smiles"]
    ].copy()

    descriptor_data["MW"] = [
        Descriptors.MolWt(mol) for mol in aryl_halide_mols
    ]

    descriptor_data["cLogP"] = [
        Descriptors.MolLogP(mol) for mol in aryl_halide_mols
    ]

    descriptor_data["TPSA"] = [
        Descriptors.TPSA(mol) for mol in aryl_halide_mols
    ]

    descriptor_data["HBA"] = [
        Descriptors.NumHAcceptors(mol) for mol in aryl_halide_mols
    ]

    descriptor_data["HBD"] = [
        Descriptors.NumHDonors(mol) for mol in aryl_halide_mols
    ]

    descriptor_data["RotatableBonds"] = [
        Descriptors.NumRotatableBonds(mol) for mol in aryl_halide_mols
    ]

    descriptor_data["FractionCSP3"] = [
        Descriptors.FractionCSP3(mol) for mol in aryl_halide_mols
    ]

    descriptor_data["HeavyAtomCount"] = [
        Descriptors.HeavyAtomCount(mol) for mol in aryl_halide_mols
    ]

    descriptor_data["RingCount"] = [
        Descriptors.RingCount(mol) for mol in aryl_halide_mols
    ]

    descriptor_data["AromaticRingCount"] = [
        Lipinski.NumAromaticRings(mol) for mol in aryl_halide_mols
    ]

    descriptor_data["MolMR"] = [
        Descriptors.MolMR(mol) for mol in aryl_halide_mols
    ]

    descriptor_data["LabuteASA"] = [
        Descriptors.LabuteASA(mol) for mol in aryl_halide_mols
    ]

    descriptor_data
    return (descriptor_data,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### PCA of molecular descriptors

    Because these descriptors have different units and numerical ranges, we first
    **standardise** them to mean 0 and standard deviation 1.

    We then use PCA to reduce the 12-dimensional descriptor space to two principal
    components for visualisation.

    The **scores** tell us where each molecule lies in this new space, while the
    **loadings** tell us which molecular descriptors contribute most strongly to
    each principal component.
    """)
    return


@app.cell(hide_code=True)
def _(descriptor_data):
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA

    descriptor_columns = [
        "MW",
        "cLogP",
        "TPSA",
        "HBA",
        "HBD",
        "RotatableBonds",
        "FractionCSP3",
        "HeavyAtomCount",
        "RingCount",
        "AromaticRingCount",
        "MolMR",
        "LabuteASA",
    ]

    X_descriptors = descriptor_data[descriptor_columns].to_numpy()

    X_descriptors_scaled = StandardScaler().fit_transform(X_descriptors)

    descriptor_pca = PCA(n_components=2)
    descriptor_scores = descriptor_pca.fit_transform(X_descriptors_scaled)
    return PCA, descriptor_columns, descriptor_pca, descriptor_scores


@app.cell(hide_code=True)
def _(descriptor_data, descriptor_pca, descriptor_scores, query_index):
    import altair as alt

    descriptor_pca_data = descriptor_data[
        ["aryl_halide_number", "aryl_halide_smiles"]
    ].copy()

    descriptor_pca_data["PC1"] = descriptor_scores[:, 0]
    descriptor_pca_data["PC2"] = descriptor_scores[:, 1]

    descriptor_pca_data["selected"] = [
        i == query_index
        for i in range(len(descriptor_pca_data))
    ]


    descriptor_pca_plot = (
        alt.Chart(descriptor_pca_data)
        .mark_circle(
            size=65,
            opacity=0.75,
        )
        .encode(
            x=alt.X(
                "PC1:Q",
                title=f"PC1 ({descriptor_pca.explained_variance_ratio_[0]*100:.1f}% variance)",
                scale=alt.Scale(zero=False),
            ),
            y=alt.Y(
                "PC2:Q",
                title=f"PC2 ({descriptor_pca.explained_variance_ratio_[1]*100:.1f}% variance)",
                scale=alt.Scale(zero=False),
            ),
            color=alt.condition(
                "datum.selected",
                alt.value("#d62728"),
                alt.value("#777777"),
            ),
            size=alt.condition(
                "datum.selected",
                alt.value(160),
                alt.value(65),
            ),
            tooltip=[
                alt.Tooltip(
                    "aryl_halide_number:Q",
                    title="Aryl halide",
                    format=".0f",
                ),
                alt.Tooltip(
                    "aryl_halide_smiles:N",
                    title="SMILES",
                ),
            ],
        )
        .properties(
            width=600,
            height=430,
            title="PCA of molecular descriptors",
        )
        .interactive()
    )

    descriptor_pca_plot
    return (alt,)


@app.cell(hide_code=True)
def _(alt, descriptor_columns, descriptor_pca, pd):
    n_loadings = 4

    loadings = pd.DataFrame(
        descriptor_pca.components_.T,
        index=descriptor_columns,
        columns=["PC1", "PC2"],
    )

    loadings

    largest_loadings = pd.concat(
        [
            loadings[pc]
            .loc[loadings[pc].abs().nlargest(n_loadings).index]
            .rename("Loading")
            .reset_index()
            .rename(columns={"index": "Descriptor"})
            .assign(PC=pc)
            for pc in ["PC1", "PC2"]
        ]
    )


    largest_loadings_plot = (
        alt.Chart(largest_loadings)
        .mark_bar()
        .encode(
            x=alt.X(
                "Loading:Q",
                title="PCA loading",
            ),
            y=alt.Y(
                "Descriptor:N",
                title=None,
                sort="-x",
            ),
            color=alt.Color(
                "PC:N",
                title=None,
                scale=alt.Scale(
                    domain=["PC1", "PC2"],
                    range=["#4c78a8", "#f58518"],
                ),
            ),
            column=alt.Column(
                "PC:N",
                title=None,
            ),
            tooltip=[
                "Descriptor:N",
                "PC:N",
                alt.Tooltip(
                    "Loading:Q",
                    format=".3f",
                ),
            ],
        )
        .properties(
            width=250,
            height=180,
            title="Largest descriptor contributions",
        )
    )

    largest_loadings_plot
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Different chemical spaces

    Both plots contain exactly the same aryl halides.

    The difference is **how the molecules were represented**:

    - **Morgan fingerprints** encode local structural environments in the molecular graph.
    - **Molecular descriptors** encode selected properties such as size, polarity,
      lipophilicity, flexibility and topology.

    Notice that molecular relationships visible in one representation may be much
    less obvious in the other.

    A molecular representation determines **which differences between molecules are
    made visible to a computational model**.

    **Chemical space is not unique: its structure depends on how we
    choose to represent the molecules.**
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
 
    """)
    return


if __name__ == "__main__":
    app.run()

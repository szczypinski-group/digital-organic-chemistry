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
        "main/workshop/data/data_table.csv"
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
        "main/workshop/ahneman-figure.png"
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
    return query_index, query_mol


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
def _(aryl_halide_mols, fp_radius, fp_size, query_index, query_mol):
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

    fp_additional_output = rdFingerprintGenerator.AdditionalOutput()
    fp_additional_output.AllocateBitInfoMap()

    display_fp = morgan_generator.GetFingerprint(
        query_mol,
        additionalOutput=fp_additional_output,
    )

    display_bit_info = fp_additional_output.GetBitInfoMap()
    active_bits = list(display_fp.GetOnBits())
    return (
        DataStructs,
        active_bits,
        aryl_halide_fps,
        display_bit_info,
        nearest_indices,
        rdFingerprintGenerator,
        similarities,
    )


@app.cell(hide_code=True)
def _(
    DataStructs,
    alt,
    aryl_halide_fps,
    fp_radius,
    fp_size,
    np,
    pd,
    query_index,
):
    fingerprint_array = np.zeros(
        fp_size.value,
        dtype=np.uint8,
    )

    DataStructs.ConvertToNumpyArray(
        aryl_halide_fps[query_index],
        fingerprint_array,
    )

    fingerprint_barcode_data = pd.DataFrame({
        "bit": np.arange(fp_size.value),
        "value": fingerprint_array,
    })

    fingerprint_barcode = (
        alt.Chart(fingerprint_barcode_data)
        .mark_rect()
        .encode(
            x=alt.X(
                "bit:O",
                title="Fingerprint bit",
                axis=None,
            ),
            color=alt.Color(
                "value:N",
                scale=alt.Scale(
                    domain=[0, 1],
                    range=["#eeeeee", "#111111"],
                ),
                legend=None,
            ),
            tooltip=[
                alt.Tooltip("bit:O", title="Bit"),
                alt.Tooltip("value:N", title="State"),
            ],
        )
        .properties(
            width=750,
            height=45,
            title=(
                f"Morgan fingerprint "
                f"(radius = {fp_radius.value}, "
                f"{fp_size.value} bits)"
            ),
        )
    )

    fingerprint_barcode
    return


@app.cell(hide_code=True)
def _(Draw, active_bits, display_bit_info, mo, query_mol):
    fingerprint_fragments_svg = Draw.DrawMorganBits(
        [
            (query_mol, bit, display_bit_info)
            for bit in active_bits
        ],
        molsPerRow=4,
        legends=[
            (
                f"Bit {bit} | r={display_bit_info[bit][0][1]}"
                + (
                    f" | {len(display_bit_info[bit])} envs"
                    if len(display_bit_info[bit]) > 1
                    else ""
                )
            )
            for bit in active_bits
        ],
        useSVG=True,
    )

    mo.Html(fingerprint_fragments_svg)
    return


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
    # Can molecular representations predict chemistry?

    So far, we have converted molecular structures into numerical representations:

    - **Morgan fingerprints**, describing structural environments;
    - **molecular descriptors**, describing selected calculated properties.

    The Buchwald–Hartwig dataset also contains an experimentally measured **reaction yield**.

    > **Can information about the molecules involved in a reaction be used to predict its yield?**

    For each reaction, we will combine representations of the:

    - aryl halide;
    - ligand;
    - base;
    - additive.

    We will then train a simple machine-learning model on some of the experimental
    reactions and ask it to predict the yields of reactions it has not seen.
    """)
    return


@app.cell(hide_code=True)
def _(Chem, DataStructs, np, pd, rdFingerprintGenerator, reactions):
    rxn_fp_size=256

    reaction_fp_generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2,
        fpSize=rxn_fp_size,
    )

    def reaction_component_fp(smiles):
        # No component present: represent with an all-zero fingerprint
        if pd.isna(smiles) or smiles == "":
            return np.zeros(rxn_fp_size, dtype=np.uint8)

        mol = Chem.MolFromSmiles(smiles)

        if mol is None:
            raise ValueError(f"Could not parse SMILES: {smiles}")

        fp = reaction_fp_generator.GetFingerprint(mol)

        fp_array = np.zeros(rxn_fp_size, dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, fp_array)

        return fp_array

    reaction_feature_rows = []

    for _, reaction in reactions.iterrows():

        has_additive = np.array(
            [0 if pd.isna(reaction["additive_smiles"]) else 1],
            dtype=np.uint8,
        )

        reaction_features = np.concatenate(
            [
                reaction_component_fp(reaction["aryl_halide_smiles"]),
                reaction_component_fp(reaction["ligand_smiles"]),
                reaction_component_fp(reaction["base_smiles"]),
                reaction_component_fp(reaction["additive_smiles"]),
                has_additive,
            ]
        )

        reaction_feature_rows.append(reaction_features)

    reaction_X = np.vstack(reaction_feature_rows)
    reaction_y = reactions["yield"].to_numpy()
    return reaction_X, reaction_y


@app.cell(hide_code=True)
def _(mo, reaction_X):
    mo.md(f"""
    ## Representing a complete reaction

    We represent each reaction by joining together the Morgan fingerprints of its
    four variable molecular components:

    \[
    X =
    FP_{{aryl\ halide}}
    \;||\;
    FP_{{ligand}}
    \;||\;
    FP_{{base}}
    \;||\;
    FP_{{additive}}
    \]

    Some experiments deliberately contain **no additive**. This is not missing
    experimental data: *no additive* is itself an experimental condition.

    For these reactions we use:

    - an all-zero additive fingerprint;
    - one additional `has_additive` feature.

    In this workshop, we will not tune the parameters of the machine learning model and will also just use a small 256-bit fingerprint.
    Each reaction is therefore represented by **{reaction_X.shape[1]} numerical
    features**, and our dataset contains **{reaction_X.shape[0]} measured reactions**.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    split_strategy = mo.ui.dropdown(
        options={
            "Random 80:20 split": "random",
            "Leave one aryl halide out": "aryl_halide",
            "Leave one ligand out": "ligand",
        },
        value="Random 80:20 split",
        label="How should we test the model?",
    )
    return (split_strategy,)


@app.cell(hide_code=True)
def _(aryl_halides, mo, reactions, split_strategy):
    aryl_holdout = mo.ui.dropdown(
        options={
            f"{int(number)}: {name}": number
            for number, name in zip(
                aryl_halides["aryl_halide_number"],
                aryl_halides["aryl_halide"],
            )
        },
        label="Aryl halide to hold out:",
        value='1: 1-chloro-4-(trifluoromethyl)benzene'
    )

    ligand_holdout = mo.ui.dropdown(
        options={
            str(name): name
            for name in sorted(reactions["ligand"].dropna().unique())
        },
        label="Ligand to hold out:",
        value=reactions["ligand"].dropna().unique()[0]
    )

    if split_strategy.value == "aryl_halide":
        holdout_selector = aryl_holdout
    elif split_strategy.value == "ligand":
        holdout_selector = ligand_holdout
    else:
        holdout_selector = mo.md(
            "The test reactions will be selected randomly."
        )

    mo.vstack(
        [
            split_strategy,
            holdout_selector,
        ]
    )
    return aryl_holdout, ligand_holdout


@app.cell(hide_code=True)
def _(
    aryl_holdout,
    ligand_holdout,
    np,
    reaction_X,
    reaction_y,
    reactions,
    split_strategy,
):
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import mean_absolute_error, r2_score

    if split_strategy.value == "random":

        (
            reaction_train_indices,
            reaction_test_indices,
        ) = train_test_split(
            np.arange(len(reactions)),
            test_size=0.20,
            random_state=42,
        )

    elif split_strategy.value == "aryl_halide":

        holdout_value = aryl_holdout.value

        reaction_test_indices = np.where(
            reactions["aryl_halide_number"].to_numpy()
            == holdout_value
        )[0]

        reaction_train_indices = np.where(
            reactions["aryl_halide_number"].to_numpy()
            != holdout_value
        )[0]

    elif split_strategy.value == "ligand":

        holdout_value = ligand_holdout.value

        reaction_test_indices = np.where(
            reactions["ligand"].to_numpy()
            == holdout_value
        )[0]

        reaction_train_indices = np.where(
            reactions["ligand"].to_numpy()
            != holdout_value
        )[0]

    reaction_X_train = reaction_X[reaction_train_indices]
    reaction_X_test = reaction_X[reaction_test_indices]

    reaction_y_train = reaction_y[reaction_train_indices]
    reaction_y_test = reaction_y[reaction_test_indices]

    yield_model = RandomForestRegressor(
        n_estimators=200,
        min_samples_leaf=2,
        n_jobs=-1,
        random_state=42,
    )

    yield_model.fit(
        reaction_X_train,
        reaction_y_train,
    )

    reaction_y_pred = yield_model.predict(
        reaction_X_test
    )

    yield_mae = mean_absolute_error(
        reaction_y_test,
        reaction_y_pred,
    )

    yield_r2 = r2_score(
        reaction_y_test,
        reaction_y_pred,
    )

    baseline_prediction = np.full_like(
        reaction_y_test,
        reaction_y_train.mean(),
        dtype=float,
    )

    baseline_mae = mean_absolute_error(
        reaction_y_test,
        baseline_prediction,
    )

    baseline_r2 = r2_score(
        reaction_y_test,
        baseline_prediction,
    )

    {
        "Random Forest MAE": round(yield_mae, 2),
        "Mean baseline MAE": round(baseline_mae, 2),
        "Random Forest R²": round(yield_r2, 3),
        "Mean baseline R²": round(baseline_r2, 3),
    }
    return reaction_y_pred, reaction_y_test


@app.cell(hide_code=True)
def _(alt, pd, reaction_y_pred, reaction_y_test):
    yield_prediction_data = pd.DataFrame(
        {
            "Experimental yield": reaction_y_test,
            "Predicted yield": reaction_y_pred,
        }
    )

    yield_range = [
        min(
            yield_prediction_data["Experimental yield"].min(),
            yield_prediction_data["Predicted yield"].min(),
        ),
        max(
            yield_prediction_data["Experimental yield"].max(),
            yield_prediction_data["Predicted yield"].max(),
        ),
    ]

    yield_points = (
        alt.Chart(yield_prediction_data)
        .mark_circle(
            size=40,
            opacity=0.55,
        )
        .encode(
            x=alt.X(
                "Experimental yield:Q",
                title="Experimental yield / %",
                scale=alt.Scale(domain=yield_range),
            ),
            y=alt.Y(
                "Predicted yield:Q",
                title="Predicted yield / %",
                scale=alt.Scale(domain=yield_range),
            ),
            tooltip=[
                alt.Tooltip("Experimental yield:Q", format=".1f"),
                alt.Tooltip("Predicted yield:Q", format=".1f"),
            ],
        )
    )

    perfect_prediction_line = (
        alt.Chart(
            pd.DataFrame(
                {
                    "Experimental yield": yield_range,
                    "Predicted yield": yield_range,
                }
            )
        )
        .mark_line(
            color="#d62728",
            strokeDash=[6, 4],
        )
        .encode(
            x="Experimental yield:Q",
            y="Predicted yield:Q",
        )
    )

    yield_prediction_plot = (
        yield_points + perfect_prediction_line
    ).properties(
        width=550,
        height=450,
        title="Predicted versus experimental reaction yield",
    ).interactive()

    yield_prediction_plot
    return


if __name__ == "__main__":
    app.run()

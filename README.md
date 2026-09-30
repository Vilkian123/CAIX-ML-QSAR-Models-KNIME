# CA IX QSAR Classification and Prospective Prediction in KNIME

This repository contains two KNIME workflows developed for ligand-based classification of **carbonic anhydrase IX (CA IX) inhibitors** using curated ChEMBL bioactivity data, RDKit molecular representations, and an **XGBoost linear classifier**.

The repository separates **model validation** from **prospective use**:

1. a reproducible **random 5-fold cross-validation workflow** used to evaluate predictive performance; and
2. a **final-model prediction workflow** that retrains the classifier on the complete curated dataset and applies it to user-supplied candidate compounds.

The prospective compound structures used in the original research project are **not included in the public repository**.

\---

## Repository contents

|File|Purpose|
|-|-|
|`XGBoostLinearClass-Random-Split-CAIX-Ki.knwf`|Data curation, molecular featurization, random 5-fold cross-validation, test-set classification metrics, and out-of-fold ROC analysis.|
|`XGBoostLinearClass-Random-Split-Test-Model-CAIX-Ki.knwf`|Retraining on the complete curated dataset and prediction of external/user-provided compounds, including class probabilities and similarity-based applicability-domain diagnostics.|



\---

## Scientific objective

CA IX is a hypoxia-associated carbonic anhydrase isoform of interest in anticancer drug discovery. The objective of these workflows is to build a reproducible classification pipeline capable of distinguishing higher-affinity from lower-affinity CA IX inhibitors from molecular structure and then applying the validated model to newly designed compounds.

The classification endpoint is derived from experimental **Ki** measurements and expressed as **pKi**:

```text
pKi = -log10(Ki \[M])
```

With Ki provided in nM, the workflow computes the equivalent transformation internally.

Compounds are classified as:

```text
Active   : median pKi >= 7.5
Inactive : median pKi < 7.5
```

A pKi threshold of 7.5 corresponds to approximately 31.6 nM.

\---

# 1\. Random 5-fold cross-validation workflow

## Workflow

`XGBoostLinearClass-Random-Split-CAIX-Ki.knwf`

This workflow is used to estimate the predictive performance of the classification model on held-out compounds.

### Main stages

```text
ChEMBL CA IX Ki data
        |
        v
Bioactivity curation
        |
        v
Duplicate aggregation / pKi calculation
        |
        v
Structural preprocessing
        |
        v
RDKit descriptors + 2048-bit fingerprint
        |
        v
Random assignment to 5 folds
        |
        v
5-fold cross-validation
        |
        +--> XGBoost Linear learner
        +--> Test-set predictor
        +--> Classification scorer
        +--> Out-of-fold probabilities
        |
        v
ROC / AUC + summary metrics
```

\---

## Bioactivity curation

The workflow starts from a ChEMBL-derived CA IX Ki table and applies a series of filters intended to produce a cleaner and more homogeneous activity dataset.

The principal curation operations include:

* retention of **exact Ki measurements** (`Standard Relation = "="`);
* retention of activities reported in **nM**;
* conversion from Ki to pKi;
* restriction to the selected ChEMBL assay format (`BAO\_0000357` in the current workflow);
* removal of compounds outside the selected AlogP range (`AlogP > 5` excluded);
* aggregation of repeated measurements by **Molecule ChEMBL ID**;
* calculation of the **median pKi** for replicated measurements;
* z-score-based pKi outlier filtering (approximately |z| <= 3);
* removal of invalid/unprocessable structures during chemical standardization.

In the current executed workflow, the final curated dataset contains **5,060 compounds** before cross-validation.

### Why median pKi?

A single compound can have multiple reported Ki measurements. Using the median reduces sensitivity to individual extreme measurements and produces one activity value per unique compound for model training.

\---

## Molecular preprocessing

Chemical structures are prepared using RDKit-based KNIME nodes. The pipeline includes:

* molecule conversion;
* salt stripping;
* structure normalization;
* canonical SMILES generation;
* RDKit molecular descriptor calculation;
* generation of a **2048-bit molecular fingerprint**;
* expansion of fingerprint bits for machine-learning input.

The feature space combines physicochemical/topological RDKit descriptors with fingerprint-derived binary features.

Examples of descriptor families represented in the workflow include:

* topological polar surface area and Labute ASA;
* hydrogen-bond donor/acceptor counts;
* rotatable bonds;
* ring and heterocycle counts;
* fraction sp3 carbon;
* Chi/Kappa topological indices;
* VSA descriptor families;
* Molecular Quantum Numbers (MQNs).

Metadata and target-related columns such as compound identifiers, pKi, fold labels, and randomization columns are excluded from model features.

> \*\*Reproducibility note:\*\* the fingerprint type, output column name, and `Expand Bit Vector` input must be kept identical across the validation and final-prediction workflows before a public release is generated.

\---

## Random 5-fold cross-validation

Random fold membership is assigned **once before the validation loop**.

The workflow uses:

```text
Random Number Assigner
        |
        | seed = 123
        v
Sorter
        v
Counter Generation
        v
Fold = mod(Counter, 5) + 1
```

This produces five approximately/equally sized folds. With 5,060 compounds, each fold contains **1,012 compounds**.

During each iteration:

```text
Fold == current fold  -> TEST
Fold != current fold  -> TRAIN
```

Therefore, every compound is used:

* once as an out-of-fold test observation;
* four times as part of a training set.

This workflow implements **random molecule-level 5-fold cross-validation**, not repeated random 80/20 holdout validation.

The folds are random but are **not explicitly stratified by activity class**.

\---

## Model

The classifier is an **XGBoost Linear Model** implemented through the KNIME XGBoost integration.

Key settings in the current workflow include:

|Parameter|Setting|
|-|-:|
|Objective|`multi:softprob`|
|Target|`Activity`|
|Boosting rounds|2000|
|L2 regularization (`lambda`)|0.5|
|L1 regularization (`alpha`)|0.0|
|Base score|0.5|
|Updater|`Shotgun`|
|Feature selector|`Cyclic`|
|Static model seed|0|

The model produces both a predicted class and class probabilities.

\---

## Classification metrics

Performance is evaluated only on the **held-out test fold** in each cross-validation iteration.

The KNIME Scorer compares:

```text
Observed class  : Activity
Predicted class : Prediction (Activity)
```

The workflow collects metrics including:

* accuracy;
* precision;
* recall / sensitivity;
* specificity;
* F-measure;
* Cohen's kappa;
* confusion-matrix counts.

Training-set metrics may also be retained as a diagnostic for overfitting, but they are **not interpreted as generalization performance**.

### Final accuracy to report

The repository should report the **mean test accuracy across the five held-out folds**.

Because each fold contains the same number of compounds in the current implementation, the mean fold accuracy is also equivalent to the accuracy computed across all pooled out-of-fold class predictions. 
> \*\*Mean 5-fold test accuracy:\*\* `0.73`  
> \*\*Cohen's kappa:\*\* `0.414`

Do not substitute the nearest-neighbour agreement score from the prospective workflow for model accuracy.

\---

## ROC curve and ROC-AUC

The active ROC node is configured with:

```text
Class column  : Activity
Positive class: 1.0
Score column  : P (Activity=1.0)
```

The ROC curve is generated from the **pooled out-of-fold predictions** produced during random 5-fold cross-validation. Thus, every probability entering the ROC calculation corresponds to a compound predicted by a model that did not use that compound for training.

This should be described as:

> \*\*ROC-AUC calculated from pooled out-of-fold predictions from random 5-fold cross-validation.\*\* 
> \*\*Pooled out-of-fold ROC-AUC:\*\* `0.794`

\---

# 2\. Final-model / prospective-prediction workflow

## Workflow

`XGBoostLinearClass-Random-Split-Test-Model-CAIX-Ki.knwf`

This workflow is intended for model application **after validation has already been performed**.

Its purpose is different from the cross-validation workflow:

```text
Complete curated dataset
        |
        v
Train final XGBoost model on all available curated compounds
        |
        +-------------------------------+
        |                               |
        |                        User candidate file
        |                               |
        |                               v
        |                    Same molecular preprocessing
        |                               |
        +---------- model --------------+
                                        |
                                        v
                                  Predictions
                                        |
                                        v
                           Tanimoto / applicability domain
```

The complete curated reference dataset is used to fit the final classifier so that no validated training information is unnecessarily discarded when generating prospective predictions.

\---

## Candidate-compound input

The original prospective structures used in the research project are intentionally **not distributed** in this repository.

The public workflow is designed to accept a user-provided compound file containing two fields:

|Column|Meaning|
|-|-|
|`Column0`|SMILES string|
|`Column1`|Compound identifier|

For clarity, users may rename these columns to `SMILES` and `Compound\_ID` in their own input file, provided the downstream reader/mapping is updated accordingly.

Example only:

```text
CCOc1ccc(...)cc1    compound\_001
CCN1CCC(...)CC1     compound\_002
```

**Do not commit unpublished or confidential candidate structures to the public repository.**

For a public GitHub/KNIME Hub version, the candidate reader should either:

* point to a small non-confidential example file; or
* remain as a user-configurable input node.

Before exporting a public `.knwf`, reset the candidate branch so that cached tables do not contain confidential SMILES or prediction results.

\---

## Final-model prediction output

The prospective workflow generates information such as:

* predicted activity class;
* probability of the inactive class;
* probability of the active class;
* maximum Tanimoto similarity to the training set;
* nearest training-set compound information;
* nearest-neighbour activity class;
* similarity-based applicability-domain label.

The exact output columns may depend on the final column-renaming configuration used in the public workflow.

\---

## Similarity-based applicability domain

For each candidate, the workflow identifies the **single most similar training compound** using Tanimoto similarity.

Current prospective-model settings include:

```text
Nearest neighbours considered: 1
Similarity metric            : Tanimoto
Applicability-domain cutoff  : >= 0.75
```

Candidates are therefore labelled approximately as:

```text
Max similarity >= 0.75 -> In-Domain
Max similarity <  0.75 -> Out-of-Domain
```

This applicability-domain definition is a pragmatic similarity-based diagnostic rather than a formal guarantee of prediction reliability.

\---

## Nearest-neighbour agreement score

The prospective workflow also contains a Scorer comparing:

```text
Prediction (Activity)
vs.
Nearest\_Train\_Activity
```

This measures **agreement between the XGBoost prediction and the activity class of the nearest training-set neighbour**.

It must **not** be interpreted as external predictive accuracy because the candidate compounds do not have experimentally measured ground-truth activities in this workflow.

The nearest-neighbour comparison is intended only as a local consistency / sanity-check diagnostic.

\---

# Data requirements

## Training data

The curation workflow expects a ChEMBL-style CA IX activity export containing the fields required by the pipeline, including relevant identifiers, chemical structures, Ki values, activity relations/units, and assay metadata.

Typical required columns include fields such as:

* Molecule ChEMBL ID;
* molecular structure / SMILES;
* Molecular Weight;
* AlogP;
* Standard Relation;
* Standard Value;
* Standard Units;
* pChEMBL Value;
* BAO Format ID;
* Target ChEMBL ID.

The exact schema should match the CSV Reader configuration in the workflow.

If the raw ChEMBL dataset is not distributed with this repository, users should download the appropriate CA IX activity data independently and configure the input node accordingly.

Third-party data remain subject to their original database terms and attribution requirements.

\---

# Software requirements

The workflows were developed in **KNIME Analytics Platform 4.7.x** and use nodes from the KNIME ecosystem including:

* KNIME base nodes;
* KNIME RDKit Integration;
* KNIME XGBoost Integration;
* KNIME JavaScript Views / ROC visualization nodes.

A newer KNIME version may automatically migrate some nodes. After migration, verify that:

1. molecular preprocessing produces the same structures;
2. the fingerprint configuration is unchanged;
3. the feature columns entering the learner are identical;
4. fold assignment still uses the fixed seed;
5. classification probabilities are mapped to the correct positive class.

\---

# How to run

## A. Cross-validation workflow

1. Import/open `XGBoostLinearClass-Random-Split-CAIX-Ki.knwf` in KNIME.
2. Configure the training-data `CSV Reader`.
3. Confirm that the reader points to the expected ChEMBL CA IX Ki dataset.
4. Execute the curation and molecular preprocessing nodes.
5. Verify that the curated dataset contains the expected number of compounds.
6. Execute the random-fold assignment component.
7. Confirm that folds 1-5 are approximately/equally populated.
8. Execute the XGBoost cross-validation component.
9. Inspect the **test Scorer**, not only the training Scorer.
10. Record mean test metrics across the five folds.
11. Execute/inspect the ROC node and record the pooled out-of-fold ROC-AUC.

## B. Final prediction workflow

1. Import/open `XGBoostLinearClass-Random-Split-Test-Model-CAIX-Ki.knwf`.
2. Configure the same curated/reference dataset used by the model pipeline.
3. Supply a private candidate file containing SMILES + compound identifiers.
4. Ensure candidate compounds undergo the same structural preprocessing and molecular representation as the training molecules.
5. Execute the final learner on the complete curated reference set.
6. Run candidate prediction.
7. Review predicted classes and probabilities.
8. Review Tanimoto similarity and applicability-domain labels before interpreting predictions.
9. Treat nearest-neighbour agreement as a diagnostic only.

\---

# Portable file paths

For public/shared workflows, avoid absolute paths such as:

```text
C:\\Users\\username\\Downloads\\...
```

Prefer KNIME's **Current workflow data area** or another relative/shared path.

A portable setup may look like:

```text
workflow-data/
├── CAIX\_training\_data.csv
└── example\_candidates.tsv
```

Generated prediction files can similarly be directed to a relative output location.

Confidential candidate files should **not** be placed in the public workflow data area.

\---

# Reproducibility

Important reproducibility settings include:

* random fold seed: **123**;
* five fixed molecule-level folds;
* XGBoost static seed enabled;
* fixed binary activity threshold (`pKi >= 7.5`);
* consistent chemical standardization between training and prospective compounds;
* identical descriptors and fingerprint representation in validation and final-model workflows.

For a reproducible publication/portfolio release, the workflows should be reset and executed from scratch after all file paths and fingerprint settings have been finalized.

\---

# Interpretation and limitations

Several limitations should be considered when interpreting model performance and prospective predictions:

1. **Random CV is not scaffold CV.**  
Random molecule-level splitting can place structurally similar compounds in both train and test folds and may therefore estimate a less demanding generalization scenario than scaffold-based validation.
2. **ChEMBL activities are heterogeneous.**  
Even after endpoint and assay-format curation, public bioactivity data can contain inter-assay and inter-laboratory variability.
3. **The activity threshold is task-specific.**  
The `pKi = 7.5` cutoff converts a continuous biochemical endpoint into a binary classification problem and should be interpreted in the context of this modeling task.
4. **Applicability-domain membership is heuristic.**  
A Tanimoto cutoff of 0.75 is a similarity-based operational definition and does not guarantee predictive correctness.
5. **Nearest-neighbour agreement is not external validation.**  
Agreement with the closest training analogue is useful context but cannot replace experimental testing.
6. **Prospective predictions require experimental validation.**  
Model output is intended to support compound prioritization, not to establish biological activity by itself.

\---



# Citation and acknowledgement

If this repository contributes to academic work, please cite the corresponding thesis/publication when available.

Bioactivity information used to construct the model was derived from **ChEMBL**, and molecular processing was performed using **RDKit** within **KNIME Analytics Platform**.

* ChEMBL: https://www.ebi.ac.uk/chembl/
* KNIME: https://www.knime.com/
* RDKit: https://www.rdkit.org/
* XGBoost: https://xgboost.ai/

\---

# License

The repository code, workflow organization, and original documentation may be distributed under the **MIT License** if that is the license selected for this repository.

Third-party software and datasets remain governed by their respective licenses and terms of use.

\---

# Author

**Maximiliano Javier Gandola Toledo** 
Drug discovery / medicinal chemistry / computational modeling

For questions about the implementation, open an issue in this repository.

\---

## Status

This repository documents a research-oriented QSAR workflow. It is intended for scientific and educational use and is **not a clinical decision-making tool**.


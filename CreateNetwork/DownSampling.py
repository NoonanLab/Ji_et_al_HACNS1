import os
import numpy as np
import scanpy as sc
import loompy as lp

BASE_DIR = "/vast/palmer/scratch/noonan/yj345/snRNA/grn"

TISSUES = ["FL", "HL", "PA1", "PA2"]

STRAT_KEY = "celltype"        
N_SUBSAMPLES = 40             
CELLS_PER_SUBSAMPLE = 3000    
RANDOM_SEED = 123             

def sample_indices(group_indices, n_total, rng):
    sizes = {ct: len(idx) for ct, idx in group_indices.items()}
    total_cells = sum(sizes.values())
    props = {ct: sizes[ct] / total_cells for ct in group_indices}

    alloc = {ct: int(np.floor(props[ct] * n_total)) for ct in group_indices}
    allocated = sum(alloc.values())
    leftover = n_total - allocated

    if leftover > 0:
        cts_sorted = sorted(group_indices.keys(), key=lambda x: props[x], reverse=True)
        for k in range(leftover):
            ct = cts_sorted[k % len(cts_sorted)]
            alloc[ct] += 1

    all_indices = []

    for ct, idx in group_indices.items():
        k = alloc[ct]

        if k == 0:
            continue

        replace = k > len(idx)
        chosen = rng.choice(idx, size=k, replace=replace)
        all_indices.append(chosen)

    all_indices = np.concatenate(all_indices)
    rng.shuffle(all_indices)

    return all_indices


def run_one_tissue(tissue):
    in_h5ad = os.path.join(
        BASE_DIR,
        tissue,
        f"adata_E10_{tissue}_filtered.h5ad"
    )

    out_dir = os.path.join(
        BASE_DIR,
        tissue,
        "sub"
    )

    os.makedirs(out_dir, exist_ok=True)

    print(f"Running tissue: {tissue}")
   
    adata = sc.read_h5ad(in_h5ad)
    
    adata.X = adata.layers["magic"].copy()
    groups = adata.obs.groupby(STRAT_KEY).groups

    group_indices = {}

    for ct, idx_labels in groups.items():
        group_indices[ct] = adata.obs.index.get_indexer(idx_labels)
  
    tissue_seed = RANDOM_SEED + TISSUES.index(tissue)
    rng = np.random.default_rng(tissue_seed)

    for i in range(N_SUBSAMPLES):
        sub_idx = sample_indices(
            group_indices,
            CELLS_PER_SUBSAMPLE,
            rng
        )

        adata_sub = adata[sub_idx, :].copy()
        
        X = adata_sub.X

        if not isinstance(X, np.ndarray):
            X = X.toarray()

        X_T = X.T

        row_attrs = {
            "Gene": np.array(adata_sub.var_names)
        }

        nGene = np.sum(X_T > 0, axis=0).flatten()
        nUMI = np.sum(X_T, axis=0).flatten()

        col_attrs = {
            "CellID": np.array(adata_sub.obs_names),
            "nGene": nGene,
            "nUMI": nUMI,
        }

        out_loom = os.path.join(
            out_dir,
            f"{tissue}_E10_magic_sub_{i + 1:02d}.loom"
        )

        lp.create(out_loom, X_T, row_attrs, col_attrs)

    print(f"\nFinished tissue: {tissue}")


def main():
    for tissue in TISSUES:
        run_one_tissue(tissue)

    print("\nAll tissues finished.")


if __name__ == "__main__":
    main()

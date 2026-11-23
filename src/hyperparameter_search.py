import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedGroupKFold
from xgboost import XGBClassifier

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import config as C

CLASSES = ["sensitivity", "neutral", "resistance"]
CLASS_INDEX = {name: position for position, name in enumerate(CLASSES)}
FOLDS = 5
SEARCH_DRAWS = [90, 91]
CANDIDATES = 40
OUT = C.DATA / "model_hyperparameter_search.csv"

SPACE = {
    "max_depth": [3, 4, 6, 8],
    "learning_rate": [0.02, 0.05, 0.1],
    "n_estimators": [200, 400, 800],
    "min_child_weight": [1, 3, 10],
    "subsample": [0.6, 0.8, 1.0],
    "colsample_bytree": [0.4, 0.6, 0.8],
    "reg_lambda": [0.5, 1.0, 5.0],
    "gamma": [0.0, 1.0],
}

DEFAULT = {
    "max_depth": 6, "learning_rate": 0.05, "n_estimators": 400, "min_child_weight": 1,
    "subsample": 0.8, "colsample_bytree": 0.8, "reg_lambda": 1.0, "gamma": 0.0,
}

X = pd.read_csv(C.DATA / "model_matrix_X.csv.gz", index_col="row_id")
y = pd.read_csv(C.DATA / "model_matrix_y.csv.gz", index_col="row_id")
FEATURES = list(X.columns)
table = X.join(y).reset_index()
print(f"{len(table):,} rows, {len(FEATURES)} features")
print(f"classes: {table.call.value_counts().to_dict()}")


def downsample(frame, draw):
    called = frame[frame.call != "neutral"]
    neutral = frame[frame.call == "neutral"]
    target = int(frame.call.value_counts().iloc[1])
    fraction = min(1.0, target / len(neutral))
    generator = np.random.default_rng(draw)
    chosen = []
    for positions in neutral.groupby(["gene", "cell_line"]).indices.values():
        wanted = min(len(positions), max(1, int(round(len(positions) * fraction))))
        chosen.append(generator.choice(positions, size=wanted, replace=False))
    return pd.concat([called, neutral.iloc[np.concatenate(chosen)]], ignore_index=True)

BLOCKS = []
for draw in SEARCH_DRAWS:
    block = downsample(table, draw)
    BLOCKS.append(
        (
            block[FEATURES].astype(float).to_numpy(),
            block.call.map(CLASS_INDEX).to_numpy(),
            block.edit.to_numpy(),
        )
    )
print(f"{len(BLOCKS)} draws prepared, {BLOCKS[0][0].shape[0]:,} rows each")


def score(parameters, seed=0):
    per_draw = []
    for features, labels, groups in BLOCKS:
        probability = np.zeros((len(labels), len(CLASSES)))
        splitter = StratifiedGroupKFold(FOLDS, shuffle=True, random_state=seed)
        
        for train, test in splitter.split(features, labels, groups=groups):
            model = XGBClassifier(
                objective="multi:softprob", num_class=len(CLASSES), tree_method="hist",
                n_jobs=-1, random_state=seed, verbosity=0, **parameters,
            )
            model.fit(features[train], labels[train])
            probability[test] = model.predict_proba(features[test])
        
        per_draw.append(
            {
                name: average_precision_score(
                    (labels == index).astype(int), probability[:, index]
                )
                for name, index in CLASS_INDEX.items()
            }
        )
    return pd.DataFrame(per_draw).mean().to_dict()


generator = np.random.default_rng(0)
sampled, seen = [DEFAULT], {tuple(sorted(DEFAULT.items()))}
while len(sampled) < CANDIDATES + 1:
    candidate = {name: generator.choice(values).item() for name, values in SPACE.items()}
    key = tuple(sorted(candidate.items()))
    if key not in seen:
        seen.add(key)
        sampled.append(candidate)

rows, started = [], time.time()

for position, parameters in enumerate(sampled):
    elapsed = time.time() - started
    result = score(parameters)
    rows.append({**parameters, **result, "is_default": parameters == DEFAULT})
    done = position + 1
    rate = (time.time() - started) / done

search = pd.DataFrame(rows)
search["macro"] = search[CLASSES].mean(axis=1)
search = search.sort_values("resistance", ascending=False)
search.to_csv(OUT, index=False)

columns = [*SPACE, *CLASSES, "macro", "is_default"]
print(search[columns].head(5).round(4).to_string(index=False))

default_row = search[search.is_default].iloc[0]
best = search.iloc[0]

print(f"\ndefault:  resistance {default_row.resistance:.4f}, macro {default_row.macro:.4f}")
print(f"best:     resistance {best.resistance:.4f}, macro {best.macro:.4f}")
print(f"gain:     {best.resistance - default_row.resistance:+.4f} resistance, "
      f"{best.macro - default_row.macro:+.4f} macro")

best_parameters = {
    name: int(best[name]) if name in {"max_depth", "n_estimators", "min_child_weight"}
    else float(best[name])
    for name in SPACE
}

print(f"\nbest parameters: {best_parameters}")

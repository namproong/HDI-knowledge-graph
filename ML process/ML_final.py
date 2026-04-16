#%%
#Data Loading & Setup
import numpy as np
import pickle
import pandas as pd

import xgboost as xgb
import lightgbm as lgb
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.neighbors import KNeighborsClassifier
 
import time
import optuna
from sklearn.model_selection import cross_val_score, train_test_split, StratifiedKFold
from sklearn.calibration import CalibratedClassifierCV
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    confusion_matrix,
    accuracy_score,
    precision_score,
    f1_score
)
#%%
from google.colab import drive
drive.mount('/content/drive')
#%%
base_path = "/content/drive/MyDrive/HDI_project/"

X_train = np.load(base_path + "X_trainHDIfinal.npy")
y_train = np.load(base_path + "y_trainHDIfinal.npy")
X_test = np.load(base_path + "X_testHDIfinal.npy")
y_test = np.load(base_path + "y_testHDIfinal.npy")
#%%
print(X_train.shape, y_train.shape) #(1460074, 1205) (1460074,)
print(X_test.shape, y_test.shape) #(365020, 1205) (365020,)
#%%
# data sampling
X_select, _, y_select, _ = train_test_split(
    X_train, y_train, train_size=0.1, stratify=y_train, random_state=42
)
#%%
# model selection 
N_ESTIMATORS = 400  
MAX_DEPTH = 8       
LEARNING_RATE = 0.1  

models = [
    # 1. Baseline: Linear Model  
    ('SGD (Logistic)', CalibratedClassifierCV(SGDClassifier(loss='log_loss', max_iter=2000, random_state=42))),
    # 2. Instance-based: KNN  
    ('KNN', KNeighborsClassifier(n_neighbors=5, n_jobs=-1)),
    # 3. Bagging: RandomForest  
    ('RandomForest', RandomForestClassifier(n_estimators=N_ESTIMATORS, max_depth=MAX_DEPTH, n_jobs=-1, random_state=42)),
    # 4. HistGradBoost  
    ('HistGradBoost', HistGradientBoostingClassifier(max_iter=N_ESTIMATORS, max_depth=MAX_DEPTH, random_state=42)),
    # 5. LightGBM  
    ('LightGBM', lgb.LGBMClassifier(
        n_estimators=N_ESTIMATORS,
        max_depth=MAX_DEPTH,
        num_leaves=2**7,  
        learning_rate=LEARNING_RATE,
        device='gpu',
        random_state=42,
        verbosity=-1
    )),
    # 6. Neural Network: MLP  
    ('MLP (Neural Net)', MLPClassifier(hidden_layer_sizes=(128, 128), max_iter=100, random_state=42)),
    # 7. XGBoost  
    ('XGBoost (Standard)', xgb.XGBClassifier(
        tree_method="hist",
        device="cuda",
        n_estimators=N_ESTIMATORS,
        max_depth=MAX_DEPTH,
        learning_rate=LEARNING_RATE,
        random_state=42
    ))
]

print(f"--- Starting Model Selection (Sample: {len(X_select)} rows) ---")
print(f"Benchmark: Estimators={N_ESTIMATORS}, Max Depth={MAX_DEPTH}, LR={LEARNING_RATE}")
#%%
results = []
for name, model in models:
    print(f"Running {name}...")
    start_time = time.time()

    # ใช้ Stratified 5-Fold CV  
    try:
        cv_scores = cross_val_score(model, X_select, y_select, cv=5, scoring='roc_auc')
        elapsed = time.time() - start_time

        results.append({
            "Model": name,
            "AUROC (Mean)": round(cv_scores.mean(), 4),
            "AUROC (Std)": round(cv_scores.std(), 4),
            "Time (sec)": round(elapsed, 2)
        })
        print(f"Done: {name} | AUC: {cv_scores.mean():.4f} | Time: {elapsed:.2f}s")
    except Exception as e:
        print(f"Error with {name}: {e}")

# สรุปผล
selection_summary = pd.DataFrame(results).sort_values(by="AUROC (Mean)", ascending=False)
print("\n=== COMPARISON SUMMARY ===")
print(selection_summary)
#%%
#Hyperparameter Tuning with Optuna
def objective(trial):

    params = {
        "n_estimators": trial.suggest_int("n_estimators", 200, 800),
        "max_depth": trial.suggest_int("max_depth", 4, 8),
        "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.2),
        "subsample": trial.suggest_float("subsample", 0.7, 1.0),
        "colsample_bytree": trial.suggest_float("colsample_bytree", 0.7, 1.0),
        "min_child_weight": trial.suggest_int("min_child_weight", 3, 10),
        "gamma": trial.suggest_float("gamma", 0, 3),
        "reg_lambda": trial.suggest_float("reg_lambda", 1, 5),

        "tree_method": "hist",
        "device": "cuda",
        "eval_metric": "auc",
        "random_state": 42,
        "early_stopping_rounds": 50
    }

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aucs = []

    for train_idx, val_idx in skf.split(X_train, y_train):

        X_tr, X_val = X_train[train_idx], X_train[val_idx]
        y_tr, y_val = y_train[train_idx], y_train[val_idx]

        model = xgb.XGBClassifier(**params)

        model.fit(
                  X_tr, y_tr,
                  eval_set=[(X_val, y_val)],
                  verbose=False
)

        pred = model.predict_proba(X_val)[:,1]
        aucs.append(roc_auc_score(y_val, pred))

    return np.mean(aucs)

#%%
study = optuna.create_study(direction="maximize", sampler=optuna.samplers.TPESampler(seed=42))
study.optimize(objective, n_trials=5)

best_params = study.best_params
print("Best params:", best_params)
#%%
best_params = {'n_estimators': 425, 'max_depth': 8, 'learning_rate': 0.14907884894416698, 'subsample': 0.8795975452591109, 'colsample_bytree': 0.7468055921327309, 'min_child_weight': 4, 'gamma': 0.17425083650459838, 'reg_lambda': 4.46470458309974}
#%%
#Out-Of-Fold (OOF) Evaluation & Threshold Optimization
kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

aucs = []
auprcs = []
oof_pred = np.zeros(len(X_train))

for fold, (train_idx, val_idx) in enumerate(kf.split(X_train, y_train)):

    X_tr, X_val = X_train[train_idx], X_train[val_idx]
    y_tr, y_val = y_train[train_idx], y_train[val_idx]

    model = xgb.XGBClassifier(
        **best_params,
        tree_method="hist",    
        device="cuda",         
        eval_metric="auc",
        random_state=42,
        early_stopping_rounds=50
    )

    model.fit(
    X_tr, y_tr,
    eval_set=[(X_val, y_val)],
    verbose=False
)

    val_prob = model.predict_proba(X_val)[:,1]

    # เก็บ OOF
    oof_pred[val_idx] = val_prob

    aucs.append(roc_auc_score(y_val, val_prob))
    auprcs.append(average_precision_score(y_val, val_prob))

# CV report
print("\n=== CV RESULT ===")
print("AUROC: {:.4f} ± {:.4f}".format(np.mean(aucs), np.std(aucs)))
print("AUPRC: {:.4f} ± {:.4f}".format(np.mean(auprcs), np.std(auprcs)))

# หา threshold จาก OOF
thresholds = np.linspace(0.1, 0.9, 200)

best_f1 = 0
best_t = 0

for t in thresholds:
    pred = (oof_pred >= t).astype(int)
    f1 = f1_score(y_train, pred)

    if f1 > best_f1:
        best_f1 = f1
        best_t = t

print("\nBest threshold (OOF):", round(best_t,3))
print("Best F1:", round(best_f1,4))

oof_auc = roc_auc_score(y_train, oof_pred)
oof_auprc = average_precision_score(y_train, oof_pred)

print("\n=== OOF RESULT ===")
print("OOF AUROC:", round(oof_auc, 4))
print("OOF AUPRC:", round(oof_auprc, 4))

#%%
#Final Model Training & Testing
final_model = xgb.XGBClassifier(
    **best_params,
    tree_method="hist",
    device="cuda",
    eval_metric="auc",
    random_state=42
)
#%%
## 1. เทรนโมเดลด้วยข้อมูล Train ทั้งหมด
final_model.fit(X_train, y_train)
print("done")
#%%
# 2. ทำนายค่า Probability ของทั้ง Train และ Test
train_prob = final_model.predict_proba(X_train)[:,1]
test_prob = final_model.predict_proba(X_test)[:,1]

print("Train AUROC:", round(roc_auc_score(y_train, train_prob), 4))
print("Train AUPRC:", round(average_precision_score(y_train, train_prob), 4))
print("Test AUROC: ", round(roc_auc_score(y_test, test_prob), 4))
print("Test AUPRC: ", round(average_precision_score(y_test, test_prob), 4))
#%%
test_pred = (test_prob >= best_t).astype(int)

cm = confusion_matrix(y_test, test_pred)
tn, fp, fn, tp = cm.ravel()

sensitivity = tp / (tp + fn)
specificity = tn / (tn + fp)
precision_val = tp / (tp + fp + 1e-10)
f1 = f1_score(y_test, test_pred)
print(f"\n=== Test Set Results (Threshold: {best_t:.3f}) ===")
print("\nConfusion Matrix:")
print(cm)
print("------------------------")
print(f"Sensitivity (Recall): {sensitivity:.4f}")
print(f"Specificity:          {specificity:.4f}")
print(f"Precision:            {precision_val:.4f}")
print(f"F1-Score:             {f1:.4f}")
#%%
#Export & Save Mode
bundle = {
    "model": final_model,
    "threshold": best_t,
    "feature_shape": X_train.shape[1]
}

with open("model_bundle.pkl", "wb") as f:
    pickle.dump(bundle, f)

print("Saved: model_bundle.pkl")

test_out = pd.DataFrame({
    "y_true": y_test,
    "probability": test_prob,
    "prediction": (test_prob >= best_t).astype(int)
})

test_out.to_csv("test_predictions.csv", index=False)
print("Saved: test_predictions.csv")
#%%
import shutil
import os

destination_path = '/content/drive/MyDrive/xgboost_output'
os.makedirs(destination_path, exist_ok=True)

shutil.move('model_bundle.pkl', os.path.join(destination_path, 'model_bundle.pkl'))
shutil.move('test_predictions.csv', os.path.join(destination_path, 'test_predictions.csv'))

print(f"Files moved to: {destination_path}")
#%%

 
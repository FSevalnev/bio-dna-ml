from sklearn.metrics import classification_report, roc_auc_score, confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import KFold
from sklearn.utils import resample
from collections import Counter
import matplotlib.pyplot as plt
import itertools
import random
import numpy as np
import os

random.seed(42)
np.random.seed(42)

def read_fasta(file):
    sequence = ""
    with open(file, "r") as gene:
        #skip the fiirst line of header
        gene.readline()
        sequence = ""
        #lower case the dna code
        for line in gene:
            sequence += line.strip().lower()
    return sequence

def clean_seq(seq):
    return "".join([c for c in seq if c in "atgc"])

def gc_count(seq):
    return (seq.count("g") + seq.count("c")) / len(seq)

def at_content(seq):
    return (seq.count("a") + seq.count("t")) / len(seq)

def seq_len(seq):
    return len(seq)

def kmer_count(seq, k=3):
    kmers = {}
    for i in range(len(seq)-k+1):
        kmer = seq[i:i+k]
        kmers[kmer] = kmers.get(kmer, 0) + 1
    return kmers

def all_kmers(k=3):
    return ["".join(p) for p in itertools.product("atgc", repeat=k)]

def kmer_features(seq, k=3):
    kmers = kmer_count(seq, k)
    all_possible = all_kmers(k)

    features = []
    for kmer in all_possible:
        features.append(kmers.get(kmer, 0) / len(seq))
    
    return features

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "..", "data")

files = {
    os.path.join(DATA_DIR, "BRCA1.txt"): 1,
    os.path.join(DATA_DIR, "TP53.txt"): 1,
    os.path.join(DATA_DIR, "EGFR.txt"): 1,
    os.path.join(DATA_DIR, "GAPDH.txt"): 0,
    os.path.join(DATA_DIR, "ACTB.txt"): 0,
}

x = []
y = []


for file, label in files.items():
    seq = clean_seq(read_fasta(file))
    
    chunk_size = 200
    
    for i in range(0, len(seq), chunk_size):
        chunk = seq[i:i+chunk_size]
        if len(chunk) == chunk_size:
            features = kmer_features(chunk, k=4)
            features.append(gc_count(chunk))
            features.append(at_content(chunk))
            x.append(features)
            y.append(label)

if len(x) == 0:
    print("ERROR: No data generated")

print("Dataset size:", len(x))
print("Feature vector length:", len(x[0]))

# разделяем классы
X_0 = [x[i] for i in range(len(y)) if y[i] == 0]
X_1 = [x[i] for i in range(len(y)) if y[i] == 1]

# downsample класса 1
X_1_down = resample(X_1, replace=False, n_samples=len(X_0), random_state=42)

X_balanced = X_0 + X_1_down
y_balanced = [0]*len(X_0) + [1]*len(X_1_down)

X_train, X_test, y_train, y_test = train_test_split(
    X_balanced, y_balanced, test_size=0.3, random_state=42, stratify=y_balanced)

model = RandomForestClassifier(
    n_estimators=100,
    random_state=42,
    class_weight="balanced"
)

model.fit(X_train, y_train)
predictions = model.predict(X_test)

importances = model.feature_importances_
kmers_list = all_kmers(4)

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
scores = cross_val_score(model, X_balanced, y_balanced, cv=cv)
pairs = sorted(zip(kmers_list, importances), key=lambda x: x[1], reverse=True)
pairs = sorted(pairs, key=lambda x: x[1], reverse=True)

print("Accuracy:", accuracy_score(y_test, predictions))

for kmer, val in pairs[:10]:
    print(kmer, val)

print(confusion_matrix(y_test, predictions))
print(classification_report(y_test, predictions))
print("CV Accuracy:", scores.mean())

probs = model.predict_proba(X_test)[:,1]
print("ROC-AUC:", roc_auc_score(y_test, probs))
print(Counter(y_balanced))

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(BASE_DIR, "..", "results")

os.makedirs(results_dir, exist_ok=True)

plt.bar([k for k, _ in pairs[:10]], [v for _, v in pairs[:10]])
plt.xticks(rotation=45)
plt.title("Top Important k-mers")
plt.savefig(os.path.join(results_dir, "kmer_plot.png"))
plt.show()

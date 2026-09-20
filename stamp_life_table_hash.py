"""Insert the survival-table fingerprint into output/canonical_results.json.

One-off for the rerun that predates the fingerprint field; results_pipeline.py now
writes it directly, so this only ever needs run once (and fails loudly if the hash
is already present and disagrees).
"""
import hashlib
import json
import os

root = os.path.dirname(os.path.abspath(__file__))
lt = os.path.join(root, "data", "life_table_surv.json")
digest = hashlib.sha256(open(lt, "rb").read()).hexdigest()[:16]
path = os.path.join(root, "output", "canonical_results.json")

doc = json.load(open(path))
assert doc.get("inputs", {}).get("life_table_sha256", digest) == digest, "hash mismatch"
doc["inputs"] = {"life_table_sha256": digest, "life_table_path": "data/life_table_surv.json"}
with open(path, "w") as fh:
    json.dump(doc, fh, indent=2)
print(f"fingerprint {digest} -> {os.path.relpath(path, root)}")

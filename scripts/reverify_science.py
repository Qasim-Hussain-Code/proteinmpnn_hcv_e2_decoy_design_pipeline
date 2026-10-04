"""Repeat verification without overwriting the frozen experiment's audit files."""
from pathlib import Path
import json
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from pipeline.common import write_tsv, write_json, sha256, read_tsv, now

def main():
    folder=ROOT/'logs/console/publication_reverification';folder.mkdir(parents=True,exist_ok=True)
    from pipeline import verification as original
    from pipeline.followup import snapshot
    snapshot()
    expected={r['path']:r['sha256'] for r in read_tsv(ROOT/'followup/original_snapshot.tsv')}
    def audited_write(path,rows,fields=None):
        # The original verification writes two audit tables. Store this repeat
        # in an ignored directory, preserving the recorded historical outputs.
        return write_tsv(folder/Path(path).name,rows,fields)
    original.write_tsv=audited_write
    start=time.monotonic(); original.verify()
    write_json(folder/'original_stage.json',dict(timestamp=now(),passed=True,seconds=time.monotonic()-start))
    from pipeline.followup_report import verify
    start=time.monotonic();verify()
    write_json(folder/'followup_stage.json',dict(timestamp=now(),passed=True,seconds=time.monotonic()-start))
    snapshot()
    assert all(sha256(ROOT/p)==digest for p,digest in expected.items())
    print('Original and follow-up scientific checks passed; original archive unchanged',flush=True)

if __name__=='__main__':main()

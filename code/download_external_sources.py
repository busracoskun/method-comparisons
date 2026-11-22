#!/usr/bin/env python3
"""Download documented external method sources into an ignored directory; no fitting."""
from pathlib import Path
import argparse,hashlib,json,urllib.request
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--download',action='store_true',help='Retrieve the pinned public source files for local use');a=p.parse_args()
 entries=json.loads((ROOT/'config/external_sources.json').read_text())
 if not a.download:
  for e in entries:print(e['method'],e['url'],e['sha256'])
  return
 for e in entries:
  data=urllib.request.urlopen(e['url'],timeout=60).read()
  if hashlib.sha256(data).hexdigest()!=e['sha256']:raise ValueError('Source checksum mismatch: '+e['file'])
  dest=ROOT/e['file'];dest.parent.mkdir(parents=True,exist_ok=True)
  if dest.exists() and dest.read_bytes()!=data:raise ValueError('Refusing to replace a different local source: '+e['file'])
  dest.write_bytes(data);print('Verified',e['file'])
if __name__=='__main__':main()

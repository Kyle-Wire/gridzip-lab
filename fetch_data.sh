#!/usr/bin/env bash
# Downloads every dataset the experiments use into data/ (git-ignored).
set -euo pipefail
cd "$(dirname "$0")" && mkdir -p data && cd data
W=https://raw.githubusercontent.com/pytorch/examples/main/word_language_model/data/wikitext-2
for f in train valid test; do [ -f $f.txt ] || curl -sSfO $W/$f.txt; done
[ -f enwik8 ] || { curl -sSfO https://mattmahoney.net/dc/enwik8.zip && unzip -q enwik8.zip && rm enwik8.zip; }
echo "docs_train.txt / docs_test.txt: the specialist set is the Markdown docs of a private"
echo "flood/heat project (29 files train, 6 held out by file); not redistributed here."

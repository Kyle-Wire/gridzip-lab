#!/bin/sh
# runs after the v9 10 MB test: v10 topics at 10 MB, then v9 on all of enwik8
cd "$(dirname "$0")"
PY=/tmp/claude-0/-home-user-climate-home-inspection/ef234dbf-5e6c-5feb-a831-f5063d80cd91/scratchpad/venv/bin/python
while pgrep -f "v8_test.py ../data/enwik8 10000000 ../results/v9_10MB" >/dev/null; do sleep 20; done
ENGINE=v10 NT=32 DECAY=0.95 TB=24 $PY v8_test.py ../data/enwik8 10000000 ../results/v10_10MB.json full > ../results/v10_10MB.log 2>&1
ENGINE=v9 TB=24 SAMPLE=5000000 $PY v8_test.py ../data/enwik8 100000000 ../results/v9_enwik8_full.json full > ../results/v9_enwik8_full.log 2>&1

#!/bin/zsh
# Keeps the QM descriptor run going to queue exhaustion, across the --hours cap and across any
# unexpected driver death. Deliberately NOT part of run_qm_descriptors.py -- the script itself is
# unmodified.
#
# OFF-SWITCH:  rm outputs/qm_descriptors/KEEP_RUNNING
#   That stops supervision only. It does NOT stop a run already in progress; to stop that, send
#   SIGTERM to the driver PID, which drains in-flight jobs and exits cleanly.
#
# Relaunches with NO --hours cap, so each launch runs until the queue is exhausted.

cd /Users/codiefreeman/OpenADMET-CYP-Blind-Challenge
PY=/Users/codiefreeman/miniconda3-arm64/envs/cyp-admet-v2/bin/python
SENTINEL=outputs/qm_descriptors/KEEP_RUNNING
SUPLOG=logs/qm_supervisor.log
MAX_RELAUNCH=20
n=0

say() { echo "$(date '+%Y-%m-%d %H:%M:%S') $1" >> $SUPLOG }
say "supervisor started (pid $$), watching for the driver to exit"

while [ -f $SENTINEL ]; do
  if ! pgrep -f "scripts/run_qm_descriptors.py" > /dev/null; then
    # driver is not running -- finished, or stopped at its cap, or died
    read -r DONE TOTAL <<< "$($PY -c "
import json
try:
    d=json.load(open('outputs/qm_descriptors/status.json')); print(d['done'], d['total_compounds'])
except Exception: print(-1, -1)")"
    if [ "$DONE" -ge "$TOTAL" ] && [ "$TOTAL" -gt 0 ]; then
      say "COMPLETE: $DONE/$TOTAL compounds. Removing sentinel and exiting."
      rm -f $SENTINEL
      break
    fi
    if [ $n -ge $MAX_RELAUNCH ]; then
      say "STOPPING: hit MAX_RELAUNCH=$MAX_RELAUNCH with $DONE/$TOTAL done. Something is wrong -- not relaunching again."
      break
    fi
    n=$((n+1))
    say "driver not running at $DONE/$TOTAL -- relaunch #$n, --jobs 6, no --hours cap"
    caffeinate -i nohup $PY scripts/run_qm_descriptors.py --jobs 6 >> logs/qm_full_run.log 2>&1 &
    sleep 90
    if pgrep -f "scripts/run_qm_descriptors.py" > /dev/null; then
      say "relaunch #$n confirmed running (pid $(pgrep -f 'scripts/run_qm_descriptors.py' | head -1))"
    else
      say "WARNING: relaunch #$n did not come up; will retry on the next poll"
    fi
  fi
  sleep 120
done
say "supervisor exiting (pid $$)"

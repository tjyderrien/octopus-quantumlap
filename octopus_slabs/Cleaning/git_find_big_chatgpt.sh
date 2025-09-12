git rev-list --objects --all | git cat-file --batch-check='%(objectname) %(objecttype) %(objectsize) %(rest)' \
  | grep blob | sort -k3 -n | tail -20

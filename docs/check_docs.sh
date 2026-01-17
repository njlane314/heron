#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

if [ ! -f "${REPO_ROOT}/Doxyfile" ]; then
    echo "Can't find ${REPO_ROOT}/Doxyfile"
    exit 1
fi

echo "Checking documentation for issues in ${REPO_ROOT}"

cd "${REPO_ROOT}"

# Run doxygen in a "strict mode" to ensure that everything is documented, this will put any issues in warnings.txt
doxygen Doxyfile &> /dev/null

# Check for issues
N_ISSUES=`wc -l docs/warnings.txt | cut -d ' ' -f 1`
if [ $N_ISSUES -eq 0 ]; then
    echo; echo "No issues found! Safe to push :)"
    exit 0
fi
echo; echo "There were $N_ISSUES issues with the documentation :("
echo "Press any key for more details..."
read -n 1

echo; echo "Raw warnings -----------------------------"
cat docs/warnings.txt

echo; echo "Interpreted warnings (may be garbled) ----"

# Process the warnings.txt file to be more human readable
COUNTER=1
while read LINE; do
    FILE=`echo $LINE | cut -d ':' -f 1`
    NUM=`echo $LINE | cut -d ':' -f 2`
    ISSUE=`echo $LINE | cut -d ':' -f 4-`

    echo; echo "Issue $COUNTER: $ISSUE"
    echo "vim +$NUM $FILE"

    COUNTER=$(($COUNTER + 1))
done < docs/warnings.txt

echo; echo "^^^ Please fix the above issues before pushing code ^^^"

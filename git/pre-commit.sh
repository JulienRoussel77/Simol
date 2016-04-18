#!/bin/sh
 
OPTIONS="--style=allman --indent=spaces=2 --convert-tabs --indent-namespaces --indent-classes --indent-col1-comments --indent-switches --indent-labels --keep-one-line-blocks --pad-oper --suffix=none"
RETURN=0
ASTYLE=$(which astyle)

if [ $? -ne 0 ]; then
    echo "[!] astyle not installed. Unable to check source file format policy." >&2
    exit 1
fi

# Redirect output to stderr.
exec 1>&2
 
# Run astyle on all .cpp/.hpp files
exec git diff --cached --name-only | grep -E "\.(cpp|hpp)" | xargs astyle $OPTIONS

exit $RETURN

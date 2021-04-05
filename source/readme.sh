#!/usr/bin/env bash

echo "[$(date +"%m-%d-%Y %H:%m:%S")] Updating README.md ..."
readme=/storage/vannostrand/software/rnaseq/README.md
sed '/.\/rnaseq -h/r'<(/storage/vannostrand/software/rnaseq/rnaseq -h) README.md > "${readme}"
echo "[$(date +"%m-%d-%Y %H:%m:%S")] Successfully Updated README.md"
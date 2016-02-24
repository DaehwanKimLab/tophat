#!/usr/bin/env bash
pbin=""
fl=$(readlink $0)
if [[ -z "$fl" ]]; then
   pbin=$(dirname $0)
 else
   pbin=$(dirname $fl)
fi
export PATH=$pbin:$PATH
$pbin/tophat "$@"

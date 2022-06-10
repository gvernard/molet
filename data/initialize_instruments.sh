#!/bin/bash

builddir=$1
search_dir=$builddir/data/instrument_data
for entry in "$search_dir"/*.json
do
  dum=`basename -s .json $entry`
  IFS='_' read -ra ADDR <<< "$dum"
  mycase=${ADDR}
  ${builddir}/bin/aux_create_instrument ${search_dir}/${mycase}_specs.json ${search_dir}/${mycase}_psf.fits #2> /dev/null
done

echo "Instruments created!"

#!/bin/sh -u

gotgnuplot=`which gnuplot 2>&1 | sed 's/.*no gnuplot in.*/X/' 2>&1 | sed 's/.*not found.*/X/'`

if [ "${gotgnuplot}" = "X" ]; then
  # Under windows, gnuplot is often named wgnuplot
  gotgnuplot=`which wgnuplot 2>&1 | sed 's/.*no wgnuplot in.*/X/' 2>&1 | sed 's/.*not found.*/X/'`
fi

if [ "${gotgnuplot}" = "X" ]; then
  echo "${0}: cannot find gnuplot in your path. Please install gnuplot"
  echo "if you wish to view graphs of the timing results."
  echo "Version 4.0 or later is required."
  exit 0
fi

perfnames=`find . -name 'time_*.res*'`
nfound=`echo ${perfnames} | wc -w`
echo "Found data for ${nfound} routines in `pwd`"
for perfname in ${perfnames}
do
  echo "Processing ${perfname}"
  ./plot_one.sh ${1-} ${perfname}
done

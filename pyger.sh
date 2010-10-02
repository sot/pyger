#!/bin/sh

SKA_bin=`dirname $0`
eval `${SKA_bin}/flt_envs -shell sh`
exec $SKA_ARCH_OS/bin/python -m pyger.__init__ "$@"

exit 1
 
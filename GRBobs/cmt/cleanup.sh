if test "${CMTROOT}" = ""; then
  CMTROOT=/data0/glast/extlib/CMT/v1r14p20031120; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
tempfile=`${CMTROOT}/mgr/cmt build temporary_name -quiet`
if test ! $? = 0 ; then tempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt -quiet cleanup -sh -pack=GRBobs -version=v0p1 -path=/data0/glast/GRBSimul/celestialSources $* >${tempfile}; . ${tempfile}
/bin/rm -f ${tempfile}


# echo "Setting GRBtemplate  in /data1/glast/GLASTRELEASE_HEAD/celestialSources"

if test "${CMTROOT}" = ""; then
  CMTROOT=/data1/glast/extlib/rh9_gcc32/CMT/v1r16p20040701; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh

tempfile=`${CMTROOT}/mgr/cmt build temporary_name -quiet`
if test ! $? = 0 ; then tempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt -quiet setup -sh -pack=GRBtemplate -version= -path=/data1/glast/GLASTRELEASE_HEAD/celestialSources  $* >${tempfile}; . ${tempfile}
/bin/rm -f ${tempfile}


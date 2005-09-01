# echo "Setting GRBobs v0p1 in /data0/glast/GRBSimul/celestialSources"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /data0/glast/extlib/CMT/v1r14p20031120
endif
source ${CMTROOT}/mgr/setup.csh

set tempfile=`${CMTROOT}/mgr/cmt build temporary_name -quiet`
if $status != 0 then
  set tempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt -quiet setup -csh -pack=GRBobs -version=v0p1 -path=/data0/glast/GRBSimul/celestialSources  $* >${tempfile}; source ${tempfile}
/bin/rm -f ${tempfile}

